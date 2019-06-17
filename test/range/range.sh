# -------------------------------------------------------------------
# Project   : Wind Field Simulation
# Command   : [runmpi=true [p=1,2,..]] [n=1,2,..] sh range.sh                                              
# Author    : michael.andre@tum.de
# Created   : 2019-06-17                                        
# ------------------------------------------------------------------- 

n=${n-10000}
p=${p-2}

if test x$runmpi != xtrue ;then
    ./range $n
else
    mpirun -np $p range $n
fi
