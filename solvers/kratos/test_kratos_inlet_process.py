import unittest
from unittest.mock import patch
from unittest.mock import MagicMock
from contextlib import AbstractContextManager
import sys

from numpy import array

sys.modules['KratosMultiphysics'] = MagicMock()
sys.modules['h5py'] = MagicMock()
import kratos_inlet_process


kratos_inlet_process.VELOCITY_X = 'VELOCITY_X'
kratos_inlet_process.VELOCITY_Y = 'VELOCITY_Y'
kratos_inlet_process.VELOCITY_Z = 'VELOCITY_Z'
kratos_inlet_process.DELTA_TIME = 'DELTA_TIME'
kratos_inlet_process.TIME = 'TIME'


class TestRegularGrid1D(unittest.TestCase):

    def setUp(self):
        self.grid = kratos_inlet_process.RegularGrid1D(1., 1., 11)
    
    def test_getitem(self):
        self.assertAlmostEqual(self.grid[0], 1.)
        self.assertAlmostEqual(self.grid[1], 1.1)
        self.assertAlmostEqual(self.grid[10], 2.)

    def test_len(self):
        self.assertEqual(len(self.grid), 11)

    def test_floor_index(self):
        self.assertEqual(self.grid.floor_index(0.999999), 0)
        self.assertEqual(self.grid.floor_index(1.15), 1)

    def test_ceil_index(self):
        self.assertEqual(self.grid.ceil_index(2.000001), 10)
        self.assertEqual(self.grid.ceil_index(1.95), 10)

    def test_has(self):
        self.assertFalse(self.grid.has(0.5))
        self.assertTrue(self.grid.has(1.))
        self.assertTrue(self.grid.has(1.5))
        self.assertTrue(self.grid.has(2.))
        self.assertFalse(self.grid.has(2.1))

    def test_lower_bound(self):
        self.assertAlmostEqual(self.grid.lower_bound, 1.)

    def test_upper_bound(self):
        self.assertAlmostEqual(self.grid.upper_bound, 2.)


class TestSubGrid1D(unittest.TestCase):

    def setUp(self):
        self.parent_grid = kratos_inlet_process.RegularGrid1D(0., 1., 11)
        self.sub_grid = kratos_inlet_process.SubGrid1D(self.parent_grid, 0.21, 0.79)

    def test_getitem(self):
        self.assertAlmostEqual(self.sub_grid[0], self.parent_grid[2])
        self.assertAlmostEqual(self.sub_grid[1], self.parent_grid[3])

    def test_len(self):
        self.assertEqual(len(self.sub_grid), 7)

    def test_has(self):
        self.assertFalse(self.sub_grid.has(0.19))
        self.assertTrue(self.sub_grid.has(0.2))
        self.assertTrue(self.sub_grid.has(0.5))
        self.assertTrue(self.sub_grid.has(0.8))
        self.assertFalse(self.sub_grid.has(0.81))

    def test_floor_index(self):
        self.assertEqual(self.sub_grid.floor_index(0.199999), 0)
        self.assertEqual(self.sub_grid.floor_index(0.35), 1)

    def test_ceil_index(self):
        self.assertEqual(self.sub_grid.ceil_index(0.800001), 6)
        self.assertEqual(self.sub_grid.ceil_index(0.65), 5)

    def test_lower_bound(self):
        self.assertAlmostEqual(self.sub_grid.lower_bound, 0.2)

    def test_upper_bound(self):
        self.assertAlmostEqual(self.sub_grid.upper_bound, 0.8)

    def test_slice(self):
        self.assertEqual(self.sub_grid.gslice, slice(2, 9))


class TestInletPanel3D(unittest.TestCase):

    def setUp(self):
        grid_x = kratos_inlet_process.RegularGrid1D(0., 2., 3)
        grid_y = kratos_inlet_process.RegularGrid1D(1., 1., 2)
        grid_z = kratos_inlet_process.RegularGrid1D(1., 1., 2)
        sub_grid_y = kratos_inlet_process.SubGrid1D(grid_y, 1., 2.)
        sub_grid_z = kratos_inlet_process.SubGrid1D(grid_z, 1., 2.)
        grid = kratos_inlet_process.Grid3D(grid_x, sub_grid_y, sub_grid_z)
        data = array(
            (
                ((0.2, 0.6), (0.1, 0.4)),
                ((0.7, 0.5), (0.2, 0.8)),
                ((0.4, 0.4), (0.5, 0.6))
            )
        )
        self.panel = kratos_inlet_process.InletPanel3D(grid, data)

    def test_update(self):
        self.panel.update(0.)
        self.assertAlmostEqual(self.panel.cut_data[0, 0], 0.2)
        self.assertAlmostEqual(self.panel.cut_data[0, 1], 0.6)
        self.assertAlmostEqual(self.panel.cut_data[1, 0], 0.1)
        self.assertAlmostEqual(self.panel.cut_data[1, 1], 0.4)
        self.panel.update(0.5)
        self.assertAlmostEqual(self.panel.cut_data[0, 0], 0.45)
        self.assertAlmostEqual(self.panel.cut_data[0, 1], 0.55)
        self.assertAlmostEqual(self.panel.cut_data[1, 0], 0.15)
        self.assertAlmostEqual(self.panel.cut_data[1, 1], 0.6)
        self.panel.update(1.)
        self.assertAlmostEqual(self.panel.cut_data[0, 0], 0.7)
        self.assertAlmostEqual(self.panel.cut_data[0, 1], 0.5)
        self.assertAlmostEqual(self.panel.cut_data[1, 0], 0.2)
        self.assertAlmostEqual(self.panel.cut_data[1, 1], 0.8)
        self.panel.update(1.5)
        self.assertAlmostEqual(self.panel.cut_data[0, 0], 0.55)
        self.assertAlmostEqual(self.panel.cut_data[0, 1], 0.45)
        self.assertAlmostEqual(self.panel.cut_data[1, 0], 0.35)
        self.assertAlmostEqual(self.panel.cut_data[1, 1], 0.7)

    def test_interpolate(self):
        self.panel.update(0.)
        self.assertAlmostEqual(self.panel.interpolate(MagicMock(Y=1., Z=1.)), 0.2)
        self.assertAlmostEqual(self.panel.interpolate(MagicMock(Y=1.5, Z=1.)), 0.15)
        self.assertAlmostEqual(self.panel.interpolate(MagicMock(Y=1., Z=1.5)), 0.4)
        self.assertAlmostEqual(self.panel.interpolate(MagicMock(Y=1.5, Z=1.5)), 0.325)


class TestInletPanel2D(unittest.TestCase):

    def setUp(self):
        grid_x = kratos_inlet_process.RegularGrid1D(0., 2., 3)
        grid_z = kratos_inlet_process.RegularGrid1D(1., 1., 2)
        sub_grid_z = kratos_inlet_process.SubGrid1D(grid_z, 1., 2.)
        grid = kratos_inlet_process.Grid2D(grid_x, sub_grid_z)
        data = array(
            (
                (0.2, 0.1),
                (0.5, 0.2),
                (0.4, 0.4)
            )
        )
        self.panel = kratos_inlet_process.InletPanel2D(grid, data)

    def test_update(self):
        self.panel.update(0.)
        self.assertAlmostEqual(self.panel.cut_data[0], 0.2)
        self.assertAlmostEqual(self.panel.cut_data[1], 0.1)
        self.panel.update(0.5)
        self.assertAlmostEqual(self.panel.cut_data[0], 0.35)
        self.assertAlmostEqual(self.panel.cut_data[1], 0.15)
        self.panel.update(1.)
        self.assertAlmostEqual(self.panel.cut_data[0], 0.5)
        self.assertAlmostEqual(self.panel.cut_data[1], 0.2)
        self.panel.update(1.5)
        self.assertAlmostEqual(self.panel.cut_data[0], 0.45)
        self.assertAlmostEqual(self.panel.cut_data[1], 0.3)

    def test_interpolate(self):
        self.panel.update(0.)
        self.assertAlmostEqual(self.panel.interpolate(MagicMock(Z=1.)), 0.2)
        self.assertAlmostEqual(self.panel.interpolate(MagicMock(Z=1.6)), 0.14)
        self.assertAlmostEqual(self.panel.interpolate(MagicMock(Z=1.9)), 0.11)


class TestLogMeanProfile(unittest.TestCase):

    def setUp(self):
        self.mean_profile = kratos_inlet_process.LogMeanProfile(0.59345, 0.02, 15.)

    def test_wind_speed(self):
        self.assertAlmostEqual(self.mean_profile.wind_speed(MagicMock(Z=0.)), 0.)
        self.assertAlmostEqual(self.mean_profile.wind_speed(MagicMock(Z=20.)), 10., places=3)

    def test_bulk_wind_speed(self):
        self.assertAlmostEqual(self.mean_profile.bulk_wind_speed, 15.)


class TestImposeWindInletProcess(unittest.TestCase):

    data2d = array(
        (
            (0.1, 0.0),
            (0.4,-0.2),
            (0.1, 0.1)
        )
    )

    data3d = array(
        (
            ((0.2,-0.1), (0.1, 0.0)),
            ((0.2, 0.0), (0.4,-0.2)),
            ((0.0, 0.1), (0.1, 0.1))
        )
    )

    class MockFile(dict, AbstractContextManager):
        pass

    @staticmethod
    def setup_settings(domain_size):
        settings = {}
        settings['inlet_model_part_name'] = 'inlet'
        settings['inlet_position'] = 0.
        settings['wind_filename'] = 'wind.h5'
        settings['ramp_time'] = 10.
        settings['z0'] = 0.
        if domain_size == 3:
            settings['y0'] = 0.
        return settings

    @staticmethod
    def setup_model(domain_size):
        if domain_size == 2:
            nodes = [MagicMock(X=0.5, Z=0.5)]
        else:
            nodes = [MagicMock(X=0.5, Y=0.5, Z=0.0)]
        model = {'inlet':
                 MagicMock(
                     Nodes=nodes,
                     ProcessInfo={
                         'DELTA_TIME': 0.1,
                         'TIME': 0.0
                     }
                 )
                }
        return model

    @staticmethod
    def setup_file(domain_size):
        mock_file = TestImposeWindInletProcess.MockFile()
        mock_file['lx'] = (2.,)
        if domain_size == 3:
            mock_file['ly'] = (1.,)
        mock_file['lz'] = (1.,)
        mock_file['log_z0'] = (0.02,)
        mock_file['z'] = (20.,)
        mock_file['umean'] = (12.3,)
        if domain_size == 2:
            mock_file['u'] = TestImposeWindInletProcess.data2d
            mock_file['w'] = TestImposeWindInletProcess.data2d
        else:
            mock_file['u'] = TestImposeWindInletProcess.data3d
            mock_file['v'] = TestImposeWindInletProcess.data3d
            mock_file['w'] = TestImposeWindInletProcess.data3d
        return mock_file

    def setup_test(self, domain_size):
        self.settings = TestImposeWindInletProcess.setup_settings(domain_size)
        self.model = TestImposeWindInletProcess.setup_model(domain_size)
        self.mock_file = TestImposeWindInletProcess.setup_file(domain_size)

    def test_init(self):
        self.setup_test(domain_size=3)
        with patch('kratos_inlet_process.ImposeWindInletProcess.OpenFile') as file_:
            file_.return_value = self.mock_file
            process = kratos_inlet_process.ImposeWindInletProcess(
                self.model, self.settings)
            self.assertEqual(process.inlet_model_part_name, 'inlet')
            self.assertEqual(process.wind_filename, 'wind.h5')
            self.assertAlmostEqual(process.inlet_position, 0.)
            self.assertAlmostEqual(process.ramp_time, 10.)
            self.assertAlmostEqual(process.lx, 2.)
            self.assertAlmostEqual(process.ly, 1.)
            self.assertAlmostEqual(process.lz, 1.)
            self.assertAlmostEqual(process.log_z0, 0.02)
            self.assertAlmostEqual(process.z, 20.)
            self.assertAlmostEqual(process.umean, 12.3)
            self.assertEqual(len(process.model_part.Nodes), 1)
            
    def test_create_log_mean_profile(self):
        self.setup_test(domain_size=3)
        with patch('kratos_inlet_process.ImposeWindInletProcess.OpenFile') as file_:
            file_.return_value = self.mock_file
            process = kratos_inlet_process.ImposeWindInletProcess(
                self.model, self.settings)
            mean_profile = process.CreateLogMeanProfile()
            node = MagicMock(Z=20.)
            self.assertAlmostEqual(mean_profile.wind_speed(node), 12.3, places=2)
    
    def test_update_inlet_position(self):
        self.setup_test(domain_size=3)
        with patch('kratos_inlet_process.ImposeWindInletProcess.OpenFile') as file_:
            file_.return_value = self.mock_file
            process = kratos_inlet_process.ImposeWindInletProcess(
                self.model, self.settings)
            process.mean_profile.bulk_wind_speed = 10.0
            process.UpdateInletPosition()
            self.assertAlmostEqual(process.inlet_position, 1.0)
            for i in range(2):
                process.UpdateInletPosition()
            self.assertAlmostEqual(process.inlet_position, 1.0)

    def test_create_3d_mappers(self):
        self.setup_test(domain_size=3)
        with patch('kratos_inlet_process.ImposeWindInletProcess.OpenFile') as file_:
            file_.return_value = self.mock_file
            process = kratos_inlet_process.ImposeWindInletProcess(
                self.model, self.settings)
            self.assertEqual(process.mappers[0][0], 'VELOCITY_X')
            self.assertEqual(process.mappers[1][0], 'VELOCITY_Y')
            self.assertEqual(process.mappers[2][0], 'VELOCITY_Z')
            mapper_u = process.mappers[0][1]
            mapper_u.update(process.inlet_position)
            self.assertAlmostEqual(mapper_u.interpolate(MagicMock(Y=0.5, Z=0.5)), 0.05)
            mapper_u.update(0.5)
            self.assertAlmostEqual(mapper_u.interpolate(MagicMock(Y=0.5, Z=0.5)), 0.075)

    def test_create_2d_mappers(self):
        self.setup_test(domain_size=2)
        with patch('kratos_inlet_process.ImposeWindInletProcess.OpenFile') as file_:
            file_.return_value = self.mock_file
            process = kratos_inlet_process.ImposeWindInletProcess(
                self.model, self.settings)
            self.assertEqual(process.mappers[0][0], 'VELOCITY_X')
            self.assertEqual(process.mappers[1][0], 'VELOCITY_Z')
            mapper_u = process.mappers[0][1]
            mapper_u.update(process.inlet_position)
            self.assertAlmostEqual(mapper_u.interpolate(MagicMock(Z=0.5)), 0.05)
            mapper_u.update(0.5)
            self.assertAlmostEqual(mapper_u.interpolate(MagicMock(Z=0.5)), 0.075)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
