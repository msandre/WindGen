# ----------------------------------------------------------------------
# Project   : Wind Field Simulation
# Command   : python wfsplot.py COMMAND [FILE]
# Author    : michael.andre@tum.de
# Created   : 2013-10-07
# ----------------------------------------------------------------------

import sys
import h5py
from pylab import *

argc = len(sys.argv)

if not argc in [2,3]:
    print("Try 'python wfsplot.py help' for more information.")
elif argc == 2:
    if sys.argv[1] == "help":
        print("Usage: python wfsplot.py COMMAND [FILE]\n",                     \
              "\n",                                                            \
              "FILE is an HDF5 formatted results file created by the\n",       \
              "Wind Field Simulation library.\n",                              \
              "\n",                                                            \
              "COMMAND      Description:\n",                                   \
              "help         print this help message and exit\n",               \
              "intensity    plot the turbulence intensity in FILE\n"           \
              "spectra      plot the spectra in FILE\n")
    else:
        print("Try 'python wfsplot.py help' for more information.")
elif argc == 3:
    cmd   = sys.argv[1]
    fname = sys.argv[2]
    fh5 = h5py.File(fname, "r")
    if cmd == "intensity":
        H = 6.0
        u = fh5.get('u')
        v = fh5.get('v')
        w = fh5.get('w')
        lz = fh5.get('lz')[0]
        z = fh5.get('z')[0]
        z0 = fh5.get('log_z0')[0]
        umean = fh5.get('umean')[0]
        utau = 0.41 * umean / log(z / z0)
        su = std(u)
        sv = std(v)
        sw = std(w)
        s = linspace(1.5*z0,lz,1000)
        uz = (utau / 0.41) * log(s / z0)
        Iu = su / uz
        Iv = sv / uz
        Iw = sw / uz
        
        fig = figure()
        ax = fig.add_subplot(111,aspect="25")
        ax.plot(Iu*100.0,s/H,"b",label=r"$I_u$")
        ax.plot(Iv*100.0,s/H,"r",label=r"$I_v$")
        ax.plot(Iw*100.0,s/H,"m",label=r"$I_w$")
        ax.set_xlim(0.0,50.0)
        ax.set_ylim(0.0,2.0)
        ax.set_xlabel(r"$\mathrm{Ia}$")
        ax.set_ylabel(r"$\mathrm{z/H}$")
        ax.xaxis.grid(color='gray', linestyle='dashed')
        ax.yaxis.grid(color='gray', linestyle='dashed')
        ax.legend(loc="upper right")
        savefig("intensity.png")
        
    elif cmd == "spectra":
        Fu_exact = lambda k1z : 52.5 * k1z / (1. + 33. * k1z)**(5./3.)
        Fv_exact = lambda k1z : 8.5 * k1z / (1. + 9.5 * k1z)**(5./3.)
        Fw_exact = lambda k1z : 1.05 * k1z / (1. + 5.3 * k1z**(5./3.))
        kxz_fit = [0.00198943678865,0.00296803190755,0.00442799361834,0.00660610400926,0.00985561722592, \
                   0.0147035515589,0.021936163255,0.0327264645157,0.0488244670342,0.0728410054812,       \
                   0.108671172505,0.162126039523,0.241875118171,0.360852414347,0.538354113993,           \
                   0.8031681112,1.19824293728,1.78765331532,2.66699204005,3.9788735773]
        Fu_fit = [0.105791883098,0.146072921739,0.197712237681,0.260742123281,0.332502406069,            \
                  0.406558248946,0.47290167654,0.520191810197,0.539347512083,0.526518073323,             \
                  0.484089036307,0.420301445225,0.347791585792,0.278950554163,0.220213310142,            \
                  0.172281158258,0.133847652859,0.103416441602,0.0795951241049,0.0611154974287]
        Fv_fit = [0.0152372457433,0.0215254521029,0.0299850234773,0.0411040986205,0.0553492277494,       \
                  0.0731107692705,0.0946668613906,0.12021012889,0.149967416411,0.184289626437,           \
                  0.223131627953,0.263635981962,0.295373762089,0.301133434932,0.274083318456,            \
                  0.227711428377,0.17965853799,0.138747111279,0.106467354064,0.0815895526812]
        Fw_fit = [0.00278550169789,0.00417130698007,0.00625080407547,0.00935650000258,0.0139334761786,   \
                  0.0205076481224,0.0295848855172,0.0414924161945,0.0562252165769,0.073363792285,        \
                  0.0920622772337,0.111017919669,0.128287788544,0.141025467993,0.146042845541,           \
                  0.141762035655,0.129420220986,0.112138914397,0.0932714460426,0.0753101851988]
        u = fh5.get('u')
        v = fh5.get('v')
        w = fh5.get('w')
        lx = fh5.get('lx')[0]
        z  = fh5.get('z')[0]
        z0 = fh5.get('z0')[0]
        umean = fh5.get('umean')[0]
        utau = 0.41 * umean / log(z / z0)
        sp = u.shape
        nx = sp[0]
        ny = sp[1]
        nz = sp[2]
        Fu = zeros(int(nx/2)+1)
        Fv = zeros(int(nx/2)+1)
        Fw = zeros(int(nx/2)+1)
        kx = array([2.0 * pi * i / lx for i in range(int(nx/2)+1)])
        kxz = kx * z / (2.0 * pi)
        count = 0
        sampling_rate = 3
        for j in range(0,ny,sampling_rate):
            for k in range(0,nz,sampling_rate):
                Fu = Fu + kx * abs(fft(u[:,j,k])[:int(nx/2)+1])**2 * lx / nx**2 / (2.0 * pi) / utau**2
                Fv = Fv + kx * abs(fft(v[:,j,k])[:int(nx/2)+1])**2 * lx / nx**2 / (2.0 * pi) / utau**2
                Fw = Fw + kx * abs(fft(w[:,j,k])[:int(nx/2)+1])**2 * lx / nx**2 / (2.0 * pi) / utau**2
                count = count + 1
        if count > 0:
            Fu = Fu / count
            Fv = Fv / count
            Fw = Fw / count
        fig = figure()
        ax = fig.add_subplot(111)
        ax.loglog(kxz, Fu, '--b', label='wfs')
        ax.loglog(kxz_fit, Fu_fit, 'or', label='Fit')
        ax.loglog(kxz, Fu_exact(kxz), 'k', label='Exact')
        ax.set_ylabel(r"$\mathrm{k_xF_u(k_x)}/u_{\tau}^2$")
        ax.set_xlabel(r"$\mathrm{fz}/U_0$")
        ax.legend(loc='lower right')
        savefig('Fu.png')
        
        fig = figure()
        ax = fig.add_subplot(111)
        ax.loglog(kxz, Fv, '--b', label='wfs')
        ax.loglog(kxz_fit, Fv_fit, 'or', label='Fit')
        ax.loglog(kxz, Fv_exact(kxz), 'k', label='Exact')
        ax.set_ylabel(r"$\mathrm{k_xF_v(k_x)}/u_{\tau}^2$")
        ax.set_xlabel(r"$\mathrm{fz}/U_0$")
        ax.legend(loc='lower right')
        savefig('Fv.png')
        
        fig = figure()
        ax = fig.add_subplot(111)
        ax.loglog(kxz, Fw, '--b', label='wfs')
        ax.loglog(kxz_fit, Fw_fit, 'or', label='Fit')
        ax.loglog(kxz, Fw_exact(kxz), 'k', label='Exact')
        ax.set_ylabel(r"$\mathrm{k_xF_w(k_x)}/u_{\tau}^2$")
        ax.set_xlabel(r"$\mathrm{fz}/U_0$")
        ax.legend(loc='lower right')
        savefig('Fw.png')
    else:
        print("COMMAND not recognized!")
        print("Try 'python wfsplot.py help' for more information.")
