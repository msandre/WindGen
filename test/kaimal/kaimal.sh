# -------------------------------------------------------------------
# Project   : Wind Field Simulation
# Command   : sh kaimal.sh 
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
conv=${conv-1}
N=${N-20}
p=${p-4}

sed -i -r -e "s/(\\{|\\s)lx\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ lx = $lx /g" kaimal.wfs
sed -i -r -e "s/(\\{|\\s)ly\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ ly = $ly /g" kaimal.wfs
sed -i -r -e "s/(\\{|\\s)lz\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ lz = $lz /g" kaimal.wfs
sed -i -r -e "s/(\\{|\\s)rx\\s*=\\s*[1-9]+[0-9]*(\\}|\\s)/ rx = $rx /g" kaimal.wfs
sed -i -r -e "s/(\\{|\\s)ry\\s*=\\s*[1-9]+[0-9]*(\\}|\\s)/ ry = $ry /g" kaimal.wfs
sed -i -r -e "s/(\\{|\\s)rz\\s*=\\s*[1-9]+[0-9]*(\\}|\\s)/ rz = $rz /g" kaimal.wfs
sed -i -r -e "s/(\\{|\\s)height\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ height = $height /g" kaimal.wfs
sed -i -r -e "s/(\\{|\\s)roughness\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ roughness = $roughness /g" kaimal.wfs
sed -i -r -e "s/(\\{|\\s)umean\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ umean = $umean /g" kaimal.wfs
sed -i -r -e "s/(\\{|\\s)conv\\s*=\\s*[1-9]+[0-9]*(\\}|\\s)/ conv = $conv /g" kaimal.wfs

if test x$donotrun != xtrue ;then
    if test x$runmpi != xtrue ;then
        for i in `seq 1 $N` ;do
            ./kaimal kaimal.wfs $i
        done
    else
        for i in `seq 1 $N` ;do
            mpirun -np $p kaimal kaimal.wfs $i
        done
    fi
fi

if cat <<EOF | python ;then :
from pylab import *
twopi = 2. * pi
z = $height
utau = ($umean) * 0.41 / log(($height) / ($roughness))
Fu_exact = lambda k1z : 52.5 * k1z / (1. + 33. * k1z)**(5./3.)
Fv_exact = lambda k1z : 8.5 * k1z / (1. + 9.5 * k1z)**(5./3.)
Fw_exact = lambda k1z : 1.05 * k1z / (1. + 5.3 * k1z**(5./3.))

kxz_fit = [
0.00198943678865,
0.00296803190755,
0.00442799361834,
0.00660610400926,
0.00985561722592,
0.0147035515589,
0.021936163255,
0.0327264645157,
0.0488244670342,
0.0728410054812,
0.108671172505,
0.162126039523,
0.241875118171,
0.360852414347,
0.538354113993,
0.8031681112,
1.19824293728,
1.78765331532,
2.66699204005,
3.9788735773
]

Fu_fit = [
0.105791883098,
0.146072921739,
0.197712237681,
0.260742123281,
0.332502406069,
0.406558248946,
0.47290167654,
0.520191810197,
0.539347512083,
0.526518073323,
0.484089036307,
0.420301445225,
0.347791585792,
0.278950554163,
0.220213310142,
0.172281158258,
0.133847652859,
0.103416441602,
0.0795951241049,
0.0611154974287
]

Fv_fit = [
0.0152372457433,
0.0215254521029,
0.0299850234773,
0.0411040986205,
0.0553492277494,
0.0731107692705,
0.0946668613906,
0.12021012889,
0.149967416411,
0.184289626437,
0.223131627953,
0.263635981962,
0.295373762089,
0.301133434932,
0.274083318456,
0.227711428377,
0.17965853799,
0.138747111279,
0.106467354064,
0.0815895526812
]

Fw_fit = [
0.00278550169789,
0.00417130698007,
0.00625080407547,
0.00935650000258,
0.0139334761786,
0.0205076481224,
0.0295848855172,
0.0414924161945,
0.0562252165769,
0.073363792285,
0.0920622772337,
0.111017919669,
0.128287788544,
0.141025467993,
0.146042845541,
0.141762035655,
0.129420220986,
0.112138914397,
0.0932714460426,
0.0753101851988
]

file = 'kaimal_0.res'
res = loadtxt(file)
n = 2**int(log2(len(res)))
x = res[:n,0]
Lx = (x[-1] - x[0])
dkx = twopi / Lx
kx = array([dkx*i for i in range(n/2+1)]) 
kxz = kx * z / twopi
Fu = zeros(n)
Fv = zeros(n)
Fw = zeros(n)

for i in range($N):
  file = 'kaimal_' + str(i) + '.res'
  res = loadtxt(file)
  res = res[:n,:]
  u = res[:,1]
  v = res[:,2]
  w = res[:,3]
  Fu = Fu + abs(fft(u))**2 * Lx / n**2 / twopi 
  Fv = Fv + abs(fft(v))**2 * Lx / n**2 / twopi 
  Fw = Fw + abs(fft(w))**2 * Lx / n**2 / twopi 

Fu = kx * Fu[:n/2+1] / ($N) / utau**2
Fv = kx * Fv[:n/2+1] / ($N) / utau**2
Fw = kx * Fw[:n/2+1] / ($N) / utau**2

ax = figure(1).add_subplot(111)
plt.hold(True)
ax.loglog(kxz, Fu, '--b', label='Simulation')
ax.loglog(kxz_fit, Fu_fit, 'or', label='Fit')
ax.loglog(kxz, Fu_exact(kxz), 'k', label='Exact')
ax.set_ylabel(r"$\mathrm{k_xF_u(k_x)}/u_{\tau}^2$")
ax.set_xlabel(r"$\mathrm{fz}/U_0$")
ax.legend(loc='lower right')
savefig('u_spectrum.png')

ax = figure(2).add_subplot(111)
plt.hold(True)
ax.loglog(kxz, Fv, '--b', label='Simulation')
ax.loglog(kxz_fit, Fv_fit, 'or', label='Fit')
ax.loglog(kxz, Fv_exact(kxz), 'k', label='Exact')
ax.set_ylabel(r"$\mathrm{k_xF_v(k_x)}/u_{\tau}^2$")
ax.set_xlabel(r"$\mathrm{fz}/U_0$")
ax.legend(loc='lower right')
savefig('v_spectrum.png')

ax = figure(3).add_subplot(111)
plt.hold(True)
ax.loglog(kxz, Fw, '--b', label='Simulation')
ax.loglog(kxz_fit, Fw_fit, 'or', label='Fit')
ax.loglog(kxz, Fw_exact(kxz), 'k', label='Exact')
ax.set_ylabel(r"$\mathrm{k_xF_w(k_x)}/u_{\tau}^2$")
ax.set_xlabel(r"$\mathrm{fz}/U_0$")
ax.legend(loc='lower right')
savefig('w_spectrum.png')
EOF
else
    exit 1
fi
