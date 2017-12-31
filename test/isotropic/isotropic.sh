# -------------------------------------------------------------------
# Project   : Wind Field Simulation
# Command   : sh isotropic.sh 
# Author    : michael.andre@tum.de
# Created   : 2013-06-21
# ------------------------------------------------------------------- 

lx=${lx-13350.}
ly=${ly-1000.}
lz=${lz-1000.}
rx=${rx-12}
ry=${ry-8}
rz=${rz-8}
height=${height-25.}
roughness=${roughness-0.037}
umean=${umean-15.89}
conv=${conv-0}
N=${N-20}
p=${p-4}

sed -i -r -e "s/(\\{|\\s)lx\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ lx = $lx /g" isotropic.wfs
sed -i -r -e "s/(\\{|\\s)ly\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ ly = $ly /g" isotropic.wfs
sed -i -r -e "s/(\\{|\\s)lz\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ lz = $lz /g" isotropic.wfs
sed -i -r -e "s/(\\{|\\s)rx\\s*=\\s*[1-9]+[0-9]*(\\}|\\s)/ rx = $rx /g" isotropic.wfs
sed -i -r -e "s/(\\{|\\s)ry\\s*=\\s*[1-9]+[0-9]*(\\}|\\s)/ ry = $ry /g" isotropic.wfs
sed -i -r -e "s/(\\{|\\s)rz\\s*=\\s*[1-9]+[0-9]*(\\}|\\s)/ rz = $rz /g" isotropic.wfs
sed -i -r -e "s/(\\{|\\s)height\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ height = $height /g" isotropic.wfs
sed -i -r -e "s/(\\{|\\s)umean\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ umean = $umean /g" isotropic.wfs
sed -i -r -e "s/(\\{|\\s)conv\\s*=\\s*[1-9]+[0-9]*(\\}|\\s)/ conv = $conv /g" isotropic.wfs

if test x$donotrun != xtrue ;then
    if test x$runmpi != xtrue ;then
        for i in `seq 1 $N` ;do
            ./isotropic isotropic.wfs $i
        done
    else
        for i in `seq 1 $N` ;do
            mpirun -np $p isotropic isotropic.wfs $i
        done
    fi
fi

if cat <<EOF | python ;then :
#using isotropic reference
from pylab import *
L = $height
u_spec_coef = (9./55.) * ($umean)**2 / ($height)**(2./3.) * L**(5./3.)
vw_spec_coef = (3./110.) * ($umean)**2 / ($height)**(2./3.) * L**(5./3.)
Fu_exact = lambda k1 : u_spec_coef / (1. + (L*k1)**2)**(5./6.)
Fvw_exact = lambda k1 : vw_spec_coef * (3. + 8.*(L*k1)**2) / (1. + (L*k1)**2)**(11./6.)

file = 'isotropic_0.res'
res = loadtxt(file)
n = 2**int(log2(len(res)))
x = res[:n,0]
Lx = (x[-1] - x[0])
dkx = 2. * pi / Lx
kx = array([dkx*i for i in range(n/2+1)]) 
Fu = zeros(n)
Fv = zeros(n)
Fw = zeros(n)

for i in range($N):
  file = 'isotropic_' + str(i) + '.res'
  res = loadtxt(file)
  res = res[:n,:]
  u = res[:,1]
  v = res[:,2]
  w = res[:,3]
  Fu = Fu + abs(fft(u))**2 * Lx / (2.*pi) / n**2
  Fv = Fv + abs(fft(v))**2 * Lx / (2.*pi) / n**2
  Fw = Fw + abs(fft(w))**2 * Lx / (2.*pi) / n**2

Fu = Fu[:n/2+1] / ($N)
Fv = Fv[:n/2+1] / ($N)
Fw = Fw[:n/2+1] / ($N)

figure(1)
plt.hold(True)
plt.loglog(kx, Fu, '--b', label='u-spectrum')
plt.loglog(kx, Fu_exact(kx), 'k', label='Exact')
savefig('u_spectrum.png')

figure(2)
plt.hold(True)
plt.loglog(kx, Fv, '--b', label='v-spectrum')
plt.loglog(kx, Fw, '-.r', label='w-spectrum')
plt.loglog(kx, Fvw_exact(kx), 'k', label='Exact')
savefig('vw_spectrum.png')

EOF
else
    exit 1
fi
