'''An inlet boundary condition process for KratosMultiphysics

license: license.txt
'''

__all__ = ['Factory', 'ImposeWindInletProcess']

from collections import namedtuple
from collections import Mapping
from math import isclose
from math import floor
from math import ceil
from math import log

import h5py

import KratosMultiphysics
from KratosMultiphysics import TIME
from KratosMultiphysics import DELTA_TIME
from KratosMultiphysics import VELOCITY_X
from KratosMultiphysics import VELOCITY_Y
from KratosMultiphysics import VELOCITY_Z
from KratosMultiphysics import Logger


class Parameters(Mapping):

    def __init__(self, kratos_parameters):
        self._kratos_parameters = kratos_parameters

    def __getitem__(self, key):
        param = self._kratos_parameters[key]
        if param.IsDouble():
            value = param.GetDouble()
        elif param.IsInt():
            value = param.GetInt()
        elif param.IsString():
            value = param.GetString()
        else:
            value = param
        return value

    def __iter__(self):
        yield from self._kratos_parameters.keys()

    def __len__(self):
        return self._kratos_parameters.size()

        
def Factory(settings, Model):
    return ImposeWindInletProcess(Model, Parameters(settings['Parameters']))


Extent = namedtuple('Extent', ['lower', 'upper'])


class RegularGrid1D:

    @property
    def lower_bound(self):
        return self.extent.lower

    @property
    def upper_bound(self):
        return self.extent.upper

    def __init__(self, start_pos, length, size):
        self.extent = Extent(start_pos, start_pos+length)
        self.size = size
        self.step_size = (self.upper_bound-self.lower_bound) / (self.size-1)

    def __getitem__(self, index):
        return self.lower_bound + self.step_size*index

    def __len__(self):
        return self.size

    def floor_index(self, coord):
        if isclose(coord, self.lower_bound, abs_tol=1e-4):
            # Guard against negative index due to floating point representation.
            return 0
        local_coord = coord - self.lower_bound
        return int(floor(local_coord / self.step_size))

    def ceil_index(self, coord):
        if isclose(coord, self.upper_bound, abs_tol=1e-4):
            return len(self) - 1
        local_coord = coord - self.lower_bound
        return int(ceil(local_coord / self.step_size))

    def has(self, coord):
        return (self.floor_index(coord) >= 0
                and self.ceil_index(coord) < len(self))


IndexSpan = namedtuple('IndexSpan', ['begin', 'end'])


class SubGrid1D:

    @property
    def lower_bound(self):
        return self.grid[self.span.begin]

    @property
    def upper_bound(self):
        return self.grid[self.span.end]

    @property
    def step_size(self):
        return self.grid.step_size

    @property
    def gslice(self):
        return slice(self.span.begin, self.span.end+1)

    def __init__(self, grid, start_pos, end_pos):
        self.grid = grid
        start_index = grid.floor_index(start_pos)
        end_index = grid.ceil_index(end_pos)
        end_index = max(end_index, start_index+1)
        self.span = IndexSpan(start_index, end_index)
        
    def __getitem__(self, index):
        return self.grid[self.span.begin + index]

    def __len__(self):
        return self.span.end - self.span.begin + 1

    def has(self, coord):
        return (self.floor_index(coord) >= 0
                and self.ceil_index(coord) < len(self))

    def floor_index(self, coord):
        if isclose(coord, self.lower_bound, abs_tol=1e-4):
            return 0
        return self.grid.floor_index(coord) - self.span.begin

    def ceil_index(self, coord):
        if isclose(coord, self.upper_bound, abs_tol=1e-4):
            return len(self) - 1
        return self.grid.ceil_index(coord) - self.span.begin


Grid3D = namedtuple('Grid3D', ['x', 'y', 'z'])


class InletPanel3D:

    def __init__(self, grid, data):
        self.grid = grid
        self.data = data
        self.dy = self.grid.y.step_size
        self.dz = self.grid.z.step_size
        self.y0 = self.grid.y.lower_bound
        self.z0 = self.grid.z.lower_bound

    def update(self, pos):
        i = self.grid.x.floor_index(pos)
        tx = (pos-self.grid.x[i]) / self.grid.x.step_size
        data_0 = self.data[i, self.grid.y.gslice, self.grid.z.gslice]
        data_1 = self.data[i+1, self.grid.y.gslice, self.grid.z.gslice]
        self.cut_data = (1.0 - tx) * data_0 + tx * data_1

    def interpolate(self, node):
        # xi and eta are local mesh cell coordinates in the interval [-1, 1].
        xi = ((node.Y - self.y0) % self.dy) / self.dy
        eta = ((node.Z - self.z0) % self.dz) / self.dz
        xi = 2.0 * (xi-0.5)
        eta = 2.0 * (eta-0.5)
        # Interpolate using bilinear shape functions.
        weights = (
            0.25 * (1.0-xi) * (1.0-eta),
            0.25 * (1.0+xi) * (1.0-eta),
            0.25 * (1.0+xi) * (1.0+eta),
            0.25 * (1.0-xi) * (1.0+eta)
        )
        j = self.grid.y.floor_index(node.Y)
        k = self.grid.z.floor_index(node.Z)
        return (
            weights[0] * self.cut_data[j, k]
            + weights[1] * self.cut_data[j+1, k]
            + weights[2] * self.cut_data[j+1, k+1]
            + weights[3] * self.cut_data[j, k+1]
        )


def weak_min(nodes, key): return key(min(nodes, key=key))


def weak_max(nodes, key): return key(max(nodes, key=key))


def get_extent(nodes, key):
    lower_bound = weak_min(nodes, key=key)
    upper_bound = weak_max(nodes, key=key)
    return Extent(lower_bound, upper_bound)


class LogMeanProfile:

    def __init__(self, friction_velocity, log_z0, bulk_wind_speed,
                 key=lambda node: node.Z):
        self.friction_velocity = friction_velocity
        self.log_z0 = log_z0
        self.bulk_wind_speed = bulk_wind_speed
        self.key = key
            
    def wind_speed(self, node):
        return (2.439 * self.friction_velocity
                * log((self.key(node) + self.log_z0) / self.log_z0))


class ImposeWindInletProcess:

    @property
    def inlet_nodes(self):
        return self.model_part.Nodes

    def __init__(self, Model, settings):
        for name, value in settings.items():
            setattr(self, name, value)
        with self.OpenFile() as file_:
            for key, value in file_.items():
                if key in ['lx', 'ly', 'lz', 'log_z0', 'z', 'umean']:
                    setattr(self, key, value[0])
        self.model_part = Model[self.inlet_model_part_name]
        self.mean_profile = self.CreateLogMeanProfile()
        if len(self.inlet_nodes) > 0:
            # In MPI we pin the file if the process has nodes on the inlet.
            # The mappers are only valid during the lifetime of the file.
            self.file_ = self.OpenFile()
            self.mappers = self.Create3DMappers(self.file_)

    def ExecuteInitialize(self):
        for node in self.inlet_nodes:
            for var, _ in self.mappers:
                node.Fix(var)
        
    def ExecuteInitializeSolutionStep(self):
        self.UpdateInletPosition()
        Logger.PrintInfo('ImposeWindInletProcess',
                         'inlet position = %e' % self.inlet_position)
        self.AssignVelocity()
        self.ApplyRamp()

    def OpenFile(self):
        return h5py.File(self.wind_filename, 'r')
        
    def CreateLogMeanProfile(self, key=lambda node: node.Z):
        z = self.z
        lz = self.lz
        log_z0 = self.log_z0
        umean = self.umean
        bulk_wind_speed = umean * (log(lz/log_z0) - 1.0) / log(z/log_z0)
        friction_velocity = umean * 0.41 / log((z + log_z0) / log_z0)
        return LogMeanProfile(friction_velocity, log_z0, bulk_wind_speed)

    def Create3DMappers(self, file_):
        nx, ny, nz = file_['u'].shape
        x_grid = RegularGrid1D(0., self.lx, nx)
        y_grid = RegularGrid1D(self.y0, self.ly, ny)
        z_grid = RegularGrid1D(self.z0, self.lz, nz)
        y_extent = get_extent(self.inlet_nodes, key=lambda node: node.Y)
        z_extent = get_extent(self.inlet_nodes, key=lambda node: node.Z)
        y_subgrid = SubGrid1D(y_grid, y_extent.lower, y_extent.upper)
        z_subgrid = SubGrid1D(z_grid, z_extent.lower, z_extent.upper)
        grid = Grid3D(x_grid, y_subgrid, z_subgrid)
        # mappers become invalid after file_ is destructed.
        mappers = ((VELOCITY_X, InletPanel3D(grid, file_['u'])),
                   (VELOCITY_Y, InletPanel3D(grid, file_['v'])),
                   (VELOCITY_Z, InletPanel3D(grid, file_['w'])))
        return mappers

    def UpdateInletPosition(self):
        dt = self.model_part.ProcessInfo[DELTA_TIME]
        self.inlet_position += dt * self.mean_profile.bulk_wind_speed
        # Use the property that the wind data is periodic to cycle the domain.
        if self.inlet_position >= self.lx:
            self.inlet_position -= self.lx

    def AssignVelocity(self):
        for var, mapper in self.mappers:
            mapper.update(self.inlet_position)
            if var == VELOCITY_X:
                for node in self.inlet_nodes:
                    vel = self.mean_profile.wind_speed(node)
                    vel += mapper.interpolate(node)
                    node.SetSolutionStepValue(var, vel)
            else:
                for node in self.inlet_nodes:
                    vel = mapper.interpolate(node)
                    node.SetSolutionStepValue(var, vel)

    def ApplyRamp(self):
        time = self.model_part.ProcessInfo[TIME]
        if time < self.ramp_time:
            scal = time / self.ramp_time
            for node in self.inlet_nodes:
                for var, _ in self.mappers:
                    vel = node.GetSolutionStepValue(var)
                    node.SetSolutionStepValue(var, scal * vel)

    def Check(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass
