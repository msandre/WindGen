# -------------------------------------------------------------------
# Project   : Wind Field Simulation
# Command   : [rx=1,2,..] [ry=1,2,..] [rz=1,2,..] [runmpi=true [p=1,2,..]] sh herm3d.sh
# Author    : michael.andre@tum.de
# Created   : 2013-06-18
# ------------------------------------------------------------------- 

rx=${rx-5}
ry=${ry-5}
rz=${rz-5}
p=${p-2}

if test x$runmpi != xtrue ;then
    ./herm3d $rx $ry $rz
else
    mpirun -np $p herm3d $rx $ry $rz
fi
