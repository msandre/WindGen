# -------------------------------------------------------------------
# Project   : Wind Field Simulation
# Command   : [rx=1,2,..] [ry=1,2,..] [runmpi=true [p=1,2,..]] sh herm2d.sh
# Author    : michael.andre@tum.de
# Created   : 2013-06-15                                        
# ------------------------------------------------------------------- 

rx=${rx-8}
ry=${ry-7}
p=${p-2}

if test x$runmpi != xtrue ;then
    ./herm2d $rx $ry
else
    mpirun -np $p herm2d $rx $ry
fi
