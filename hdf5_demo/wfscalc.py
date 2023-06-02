# ----------------------------------------------------------------------
# Project   : Wind Field Simulation
# Command   : python wfscalc.py COMMAND [FILE]
# Author    : michael.andre@tum.de
# Created   : 2013-09-30
# ----------------------------------------------------------------------

import sys
import h5py
from numpy import *

argc = len(sys.argv)

if not argc in [2,3]:
    print("Try 'python wfscalc.py help' for more information.")
elif argc == 2:
    if sys.argv[1] == "help":
        print("Usage: python wfscalc.py COMMAND [FILE]\n",                     \
              "\n",                                                            \
              "FILE is an HDF5 formatted results file created by the\n",       \
              "Wind Field Simulation library.\n",                              \
              "\n",                                                            \
              "COMMAND      Description:\n",                                   \
              "help         print this help message and exit\n",               \
              "spectra      print the name of the spectra model in FILE\n",    \
              "dimension    print the domain dimensions in FILE\n",            \
              "z0           print the roughness length in FILE\n",             \
              "logz0        print the log roughness length in FILE\n",         \
              "z            print the reference height in FILE\n",             \
              "velmax       print the max velocity fluctuations in FILE\n",    \
              "velmin       print the min velocity fluctuations in FILE\n",    \
              "umean        print the mean velocity U(z) in FILE\n",           \
              "ubulk        print the mean bulk velocity in FILE\n",           \
              "stddev       print the standard deviations in FILE\n"           \
              "correlation  print the correlation coefficients in FILE\n")    
    else:
        print("Try 'python wfscalc.py help' for more information.")
elif argc == 3:
    cmd   = sys.argv[1]
    fname = sys.argv[2]
    fh5 = h5py.File(fname, "r")
    if cmd == "spectra":
        print(fh5.get('spectra')[0])
    elif cmd == "dimension":
        print("lx = ", fh5.get('lx')[0])
        print("ly = ", fh5.get('ly')[0])
        print("lz = ", fh5.get('lz')[0])
    elif cmd == "z0":
        print("z0 = ", fh5.get('z0')[0])
    elif cmd == "logz0":
        print("logz0 = ", fh5.get('log_z0')[0])
    elif cmd == "z":
        print("z = ", fh5.get('z')[0])
    elif cmd == "velmax":
        u = fh5.get('u')
        v = fh5.get('v')
        w = fh5.get('w')
        print("u-fluctuation maximum = ", amax(u))
        print("v-fluctuation maximum = ", amax(v))
        print("w-fluctuation maximum = ", amax(w))
    elif cmd == "velmin":
        u = fh5.get('u')
        v = fh5.get('v')
        w = fh5.get('w')
        print("u-fluctuation minimum = ", amin(u))
        print("v-fluctuation minimum = ", amin(v))
        print("w-fluctuation minimum = ", amin(w))
    elif cmd == "umean":
        print("umean = ", fh5.get('umean')[0])
    elif cmd == "ubulk":
        umean = fh5.get('umean')[0]
        z = fh5.get('z')[0]
        z0 = fh5.get('log_z0')[0]
        lz = fh5.get('lz')[0]
        ubulk = umean * (log(lz/z0) - 1.0) / log(z/z0)
        print("ubulk = ", ubulk)
    elif cmd == "stddev":
        varu = []
        varv = []
        varw = []
        ny = fh5.get('u').shape[1]
        nz = fh5.get('u').shape[2]
        for j in range(0,ny,ny/8):
            for k in range(0,nz,nz/8):
                varu.append(var(fh5.get('u')[:,j,k]))
                varv.append(var(fh5.get('v')[:,j,k]))
                varw.append(var(fh5.get('w')[:,j,k]))
        print("u-standard deviation = ", sqrt(mean(varu)))
        print("v-standard deviation = ", sqrt(mean(varv))) 
        print("w-standard deviation = ", sqrt(mean(varw)))
        print("v-variance / u-variance = ", mean(varv) / mean(varu))
        print("w-variance / u-variance = ", mean(varw) / mean(varu))
    elif cmd == "correlation":
        sp = fh5.get('u').shape
        n = sp[0] * sp[1] * sp[2]
        u = fh5.get('u')
        v = fh5.get('v')
        w = fh5.get('w')
        m = corrcoef([u[:,:,:].reshape((n)), v[:,:,:].reshape((n)), w[:,:,:].reshape((n))])
        print("uv-correlation coefficient = ", m[0,1])
        print("uw-correlation coefficient = ", m[0,2])
        print("vw-correlation coefficient = ", m[1,2])
    else:
        print("COMMAND not recognized!")
        print("Try 'python wfscalc.py help' for more information.")
