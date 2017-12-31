# -------------------------------------------------------------------
# Project   : Wind Field Simulation
# Command   : sh hyp2f1.sh 
# Author    : michael.andre@tum.de
# Created   : 2013-06-24
# ------------------------------------------------------------------- 

if cat <<EOF | python ;then :
import subprocess
from numpy import *
from scipy import special

x = linspace(-10.,1e-2, 100)
for z in x:
    val = float(subprocess.check_output(["./hyp2f1", "{0:17.14f}".format(z)]))
    err = abs(special.hyp2f1(1./3., 17./6., 4./3., z) - val)
    print "hyp2f1({0:f}) = {1:f},  Error = {2:e}".format(z, val, err)
EOF
else
    exit 1
fi
