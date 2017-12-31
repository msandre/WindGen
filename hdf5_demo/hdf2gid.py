########################################################################
#
# Make a GiD post file for visualizing velocity fluctuations.
#
########################################################################

import h5py
from pylab import *
import sys

name = sys.argv[1]

fh5 = h5py.File(name, 'r')
lx  = fh5.get('lx')[0]
ly  = fh5.get('ly')[0]
lz  = fh5.get('lz')[0]
sp  = fh5.get('u').shape
nx  = sp[0]
ny  = sp[1]
nz  = sp[2]
ix  = nx/64
jx  = ny
kx  = nz
u   = fh5.get('u')[:ix,:jx,:kx]
v   = fh5.get('v')[:ix,:jx,:kx]
w   = fh5.get('w')[:ix,:jx,:kx]
nid = zeros((ix,jx,kx),dtype=int)
dx = lx / float(nx - 1)
dy = ly / float(ny - 1)
dz = lz / float(nz - 1)

name = name.split('.')[0]
f   = open(name + '.post.msh', 'w')

f.write('MESH "wfs" dimension 3 Elemtype Hexahedra Nnode 8\n')
f.write('Coordinates\n')
f.write('# node number coordinate_x coordinate_y coordinate_z\n')

num = 1
x = 0.0
for i in range(ix):
    y = 0.0
    for j in range(jx):
        z = 0.0
        for k in range(kx):
            f.write('{0:d} {1:f} {2:f} {3:f}\n'.format(num, x, y, z))
            nid[i,j,k] = num
            num = num + 1
            z = z + dz
        y = y + dy
    x = x + dx

f.write('End Coordinates\n')
f.write('Elements\n')
f.write('# element node_1 node_2 node_3 node_4 node_5 node_6 node_7 node_8\n')
num = 1
for i in range(ix-1):
    for j in range(jx-1):
        for k in range(kx-1):
            n1 = nid[i,j,k]
            n2 = nid[i+1,j,k]
            n3 = nid[i+1,j+1,k]
            n4 = nid[i,j+1,k]
            n5 = nid[i,j,k+1]
            n6 = nid[i+1,j,k+1]
            n7 = nid[i+1,j+1,k+1]
            n8 = nid[i,j+1,k+1]
            f.write('{0:d} {1:d} {2:d} {3:d} {4:d} {5:d} {6:d} {7:d} {8:d}\n'.format(num,n1,n2,n3,n4,n5,n6,n7,n8))
            num = num + 1

f.write('End Elements\n')
f.close()

f = open(name + '.post.res', 'w')
f.write('GiD Post Results File 1.0\n')
f.write('Result "Fluctuations" "Wind Analysis" 1 Vector OnNodes\n')
f.write('ComponentNames "u", "v", "w"\n')
f.write('Values\n')

num = 1
for i in range(ix):
    for j in range(jx):
        for k in range(kx):
            f.write('{0:d} {1:f} {2:f} {3:f}\n'.format(num, u[i,j,k], v[i,j,k], w[i,j,k]))
            num = num + 1

f.write('End Values\n')
f.close()
fh5.close()
