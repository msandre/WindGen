# -------------------------------------------------------------------
# Project   : Wind Field Simulation
# Command   : [runmpi=true [p=1,2,..]] [n=1,2,..] sh rv.sh                                              
# Author    : michael.andre@tum.de
# Created   : 2013-06-12                                        
# ------------------------------------------------------------------- 

n=${n-50000}
p=${p-2}

if test x$runmpi != xtrue ;then
    ./rv $n
else
    mpirun -np $p rv $n
fi

if cat <<EOF | python ;then :
from pylab import *
rvs = loadtxt('rv.res')

figure(1)
H, xedges, yedges = histogram2d(rvs[:,0], rvs[:,1], bins=50, range=[[-3.,3.],[-3.,3.]], normed=1)
x, y = np.meshgrid(xedges, yedges)
z = e**(-x**2 - y**2) / pi 
extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
plt.hold(True)
plt.imshow(H, extent=extent, interpolation='nearest')
plt.colorbar()
cs = plt.contour(x, y, z, colors='k')
plt.clabel(cs, inline=1, fontsize=10)
savefig('rv.png')
EOF
else
    exit 1
fi
