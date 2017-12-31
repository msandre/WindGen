# -------------------------------------------------------------------
# Project   : Wind Field Simulation
# Command   : sh eigen.sh 
# Author    : michael.andre@tum.de
# Created   : 2013-08-20
# ------------------------------------------------------------------- 

if cat <<EOF | python ;then :
import subprocess
from numpy import *
from scipy import special

k1 = random.normal() * 100.
k2 = random.normal() * 100.
k3 = random.normal() * 100.
G  = 3.8
L  = 15.
LL = L**2

def phi11(k1, k2, k3):
    kk1 = k1 * k1
    kk2 = k2 * k2
    kk3 = k3 * k3
    kk  = kk1 + kk2 + kk3
    beta = G * (kk*LL)**(-0.333333) / sqrt(special.hyp2f1(1./3., 17./6., 4./3., -1. / (kk*LL)))
    k30 = k3  + beta * k1
    kk0 = kk1 + kk2 + k30 * k30
    top = beta * k1 * sqrt(kk1 + kk2)
    bot = (kk0 - k30*k1*beta)
    c1 = beta * kk1 * ( kk0 - 2.*k30*k30 + beta*k1*k30 ) / ( kk * (kk1 + kk2) )
    c2 = k2 * kk0 * arctan2(top,bot) / (kk1 + kk2)**1.5
    zeta1 = c1 - k2 * c2 / k1
    zeta2 = k2 * c1 / k1 + c2
    return kk0 - kk1 - 2.*k1*k30*zeta1 + (kk1 + kk2)*zeta1**2

def phi22(k1, k2, k3):
    kk1 = k1 * k1
    kk2 = k2 * k2
    kk3 = k3 * k3
    kk  = kk1 + kk2 + kk3
    beta = G * (kk*LL)**(-0.333333) / sqrt(special.hyp2f1(1./3., 17./6., 4./3., -1. / (kk*LL)))
    k30 = k3  + beta * k1
    kk0 = kk1 + kk2 + k30 * k30
    top = beta * k1 * sqrt(kk1 + kk2)
    bot = (kk0 - k30*k1*beta)
    c1 = beta * kk1 * ( kk0 - 2.*k30*k30 + beta*k1*k30 ) / ( kk * (kk1 + kk2) )
    c2 = k2 * kk0 * arctan2(top,bot) / (kk1 + kk2)**1.5
    zeta1 = c1 - k2 * c2 / k1
    zeta2 = k2 * c1 / k1 + c2
    return kk0 - kk2 - 2.*k2*k30*zeta2 + (kk1 + kk2)*zeta2**2

def phi33(k1, k2, k3):
    kk1 = k1 * k1
    kk2 = k2 * k2
    kk3 = k3 * k3
    kk  = kk1 + kk2 + kk3
    beta = G * (kk*LL)**(-0.333333) / sqrt(special.hyp2f1(1./3., 17./6., 4./3., -1. / (kk*LL)))
    k30 = k3  + beta * k1
    kk0 = kk1 + kk2 + k30 * k30
    top = beta * k1 * sqrt(kk1 + kk2)
    bot = (kk0 - k30*k1*beta)
    c1 = beta * kk1 * ( kk0 - 2.*k30*k30 + beta*k1*k30 ) / ( kk * (kk1 + kk2) )
    c2 = k2 * kk0 * arctan2(top,bot) / (kk1 + kk2)**1.5
    zeta1 = c1 - k2 * c2 / k1
    zeta2 = k2 * c1 / k1 + c2
    return (kk1 + kk2) * kk0**2 / kk**2

def phi12(k1, k2, k3):
    kk1 = k1 * k1
    kk2 = k2 * k2
    kk3 = k3 * k3
    kk  = kk1 + kk2 + kk3
    beta = G * (kk*LL)**(-0.333333) / sqrt(special.hyp2f1(1./3., 17./6., 4./3., -1. / (kk*LL)))
    k30 = k3  + beta * k1
    kk0 = kk1 + kk2 + k30 * k30
    top = beta * k1 * sqrt(kk1 + kk2)
    bot = (kk0 - k30*k1*beta)
    c1 = beta * kk1 * ( kk0 - 2.*k30*k30 + beta*k1*k30 ) / ( kk * (kk1 + kk2) )
    c2 = k2 * kk0 * arctan2(top,bot) / (kk1 + kk2)**1.5
    zeta1 = c1 - k2 * c2 / k1
    zeta2 = k2 * c1 / k1 + c2
    return -k1*k2  - k1*k30*zeta2 - k2*k30*zeta1 + (kk1 + kk2)*zeta1*zeta2

def phi13(k1, k2, k3):
    kk1 = k1 * k1
    kk2 = k2 * k2
    kk3 = k3 * k3
    kk  = kk1 + kk2 + kk3
    beta = G * (kk*LL)**(-0.333333) / sqrt(special.hyp2f1(1./3., 17./6., 4./3., -1. / (kk*LL)))
    k30 = k3  + beta * k1
    kk0 = kk1 + kk2 + k30 * k30
    top = beta * k1 * sqrt(kk1 + kk2)
    bot = (kk0 - k30*k1*beta)
    c1 = beta * kk1 * ( kk0 - 2.*k30*k30 + beta*k1*k30 ) / ( kk * (kk1 + kk2) )
    c2 = k2 * kk0 * arctan2(top,bot) / (kk1 + kk2)**1.5
    zeta1 = c1 - k2 * c2 / k1
    zeta2 = k2 * c1 / k1 + c2
    return -k1*k30 + (kk1 + kk2)*zeta1

def phi23(k1, k2, k3):
    kk1 = k1 * k1
    kk2 = k2 * k2
    kk3 = k3 * k3
    kk  = kk1 + kk2 + kk3
    beta = G * (kk*LL)**(-0.333333) / sqrt(special.hyp2f1(1./3., 17./6., 4./3., -1. / (kk*LL)))
    k30 = k3  + beta * k1
    kk0 = kk1 + kk2 + k30 * k30
    top = beta * k1 * sqrt(kk1 + kk2)
    bot = (kk0 - k30*k1*beta)
    c1 = beta * kk1 * ( kk0 - 2.*k30*k30 + beta*k1*k30 ) / ( kk * (kk1 + kk2) )
    c2 = k2 * kk0 * arctan2(top,bot) / (kk1 + kk2)**1.5
    zeta1 = c1 - k2 * c2 / k1
    zeta2 = k2 * c1 / k1 + c2
    return -k2*k30 + (kk1 + kk2)*zeta2

a = matrix([[phi11(k1,k2,k3), phi12(k1,k2,k3), phi13(k1,k2,k3)], \
            [phi12(k1,k2,k3), phi22(k1,k2,k3), phi23(k1,k2,k3)], \
            [phi13(k1,k2,k3), phi23(k1,k2,k3), phi33(k1,k2,k3)]])

print a
w, v = linalg.eig(a)
print 'Numpy eigenvalues = ({0:f}, {1:f}, {2:f})'.format(w[0], w[1], w[2])
print 'Numpy eigenvector = ({0:f}, {1:f}, {2:f})'.format(v[0,0], v[1,0], v[2,0])
print 'Numpy eigenvector = ({0:f}, {1:f}, {2:f})'.format(v[0,1], v[1,1], v[2,1])
print 'Numpy eigenvector = ({0:f}, {1:f}, {2:f})'.format(v[0,2], v[1,2], v[2,2])

subprocess.call(["./eigen", "{0:17.10e}".format(a[0,0]), "{0:17.10e}".format(a[0,1]), "{0:17.10e}".format(a[0,2]), \
                 "{0:17.10e}".format(a[1,0]), "{0:17.10e}".format(a[1,1]), "{0:17.10e}".format(a[1,2]), \
                 "{0:17.10e}".format(a[2,0]), "{0:17.10e}".format(a[2,1]), "{0:17.10e}".format(a[2,2])])
EOF
else
    exit 1
fi

