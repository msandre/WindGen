# -------------------------------------------------------------------
# Project   : Wind Field Simulation
# Command   : sh simiu.sh 
# Author    : michael.andre@tum.de
# Created   : 2013-06-24
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

sed -i -r -e "s/(\\{|\\s)lx\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ lx = $lx /g" simiu.wfs
sed -i -r -e "s/(\\{|\\s)ly\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ ly = $ly /g" simiu.wfs
sed -i -r -e "s/(\\{|\\s)lz\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ lz = $lz /g" simiu.wfs
sed -i -r -e "s/(\\{|\\s)rx\\s*=\\s*[1-9]+[0-9]*(\\}|\\s)/ rx = $rx /g" simiu.wfs
sed -i -r -e "s/(\\{|\\s)ry\\s*=\\s*[1-9]+[0-9]*(\\}|\\s)/ ry = $ry /g" simiu.wfs
sed -i -r -e "s/(\\{|\\s)rz\\s*=\\s*[1-9]+[0-9]*(\\}|\\s)/ rz = $rz /g" simiu.wfs
sed -i -r -e "s/(\\{|\\s)height\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ height = $height /g" simiu.wfs
sed -i -r -e "s/(\\{|\\s)roughness\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ roughness = $roughness /g" simiu.wfs
sed -i -r -e "s/(\\{|\\s)umean\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)/ umean = $umean /g" simiu.wfs
sed -i -r -e "s/(\\{|\\s)conv\\s*=\\s*[1-9]+[0-9]*(\\}|\\s)/ conv = $conv /g" simiu.wfs

if test x$donotrun != xtrue ;then
    if test x$runmpi != xtrue ;then
	for i in `seq 1 $N` ;do
	    ./simiu simiu.wfs $i
	done
    else
	for i in `seq 1 $N` ;do
	    mpirun -np $p simiu simiu.wfs $i
	done
    fi
fi

if cat <<EOF | python ;then :
from pylab import *
twopi = 2. * pi
z = $height
utau = ($umean) * 0.41 / log(($height) / ($roughness))
Fu_exact = lambda k1z : 100. * k1z / (1. + 50. * k1z)**(5./3.)
Fv_exact = lambda k1z : 7.5 * k1z / (1. + 9.5 * k1z)**(5./3.)
Fw_exact = lambda k1z : 1.68 * k1z / (1. + 10. * k1z**(5./3.))

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
0.135076768071,
0.184170322046,
0.245212735778,
0.316548781923,
0.392916298653,
0.465079505735,
0.521427574086,
0.551413812678,
0.549089690412,
0.51487548801,
0.455537405475,
0.383096877338,
0.310975214435,
0.247606556782,
0.194924802707,
0.152083678101,
0.117810912628,
0.09080428043,
0.0697738210809,
0.0535239829467
]

Fv_fit = [
0.0202630000695,
0.0283613220336,
0.0390895071421,
0.0529528235206,
0.070396451349,
0.0917589211407,
0.117276373664,
0.147182300596,
0.18184796016,
0.221531045986,
0.2645517455,
0.302642618769,
0.318454665522,
0.299413117773,
0.254585565134,
0.203382427932,
0.157931537051,
0.121447650363,
0.0931467118796,
0.0713919673776
]

Fw_fit = [
0.00397033357727,
0.00594901310996,
0.00891107136335,
0.0133013005026,
0.0196715692813,
0.0285973765686,
0.0405203740584,
0.0555685571604,
0.0734341312086,
0.0933379794721,
0.114011368126,
0.133549030797,
0.149084786316,
0.157018535199,
0.154908433633,
0.1433612154,
0.125518753343,
0.10517906864,
0.0853503854154,
0.0677578409336
]

file = 'simiu_0.res'
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
  file = 'simiu_' + str(i) + '.res'
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
ax.loglog(kxz, Fu, '--b', label='u-spectrum')
ax.loglog(kxz_fit, Fu_fit, 'or', label='Fit')
ax.loglog(kxz, Fu_exact(kxz), 'k', label='Exact')
ax.set_ylabel(r"$\mathrm{k_xF_u(k_x)}/u_{\tau}^2$")
ax.set_xlabel(r"$\mathrm{fz}/U_0$")
savefig('u_spectrum.png')

ax = figure(2).add_subplot(111)
plt.hold(True)
ax.loglog(kxz, Fv, '--b', label='v-spectrum')
ax.loglog(kxz_fit, Fv_fit, 'or', label='Fit')
ax.loglog(kxz, Fv_exact(kxz), 'k', label='Exact')
ax.set_ylabel(r"$\mathrm{k_xF_v(k_x)}/u_{\tau}^2$")
ax.set_xlabel(r"$\mathrm{fz}/U_0$")
savefig('v_spectrum.png')

ax = figure(3).add_subplot(111)
plt.hold(True)
ax.loglog(kxz, Fw, '--b', label='w-spectrum')
ax.loglog(kxz_fit, Fw_fit, 'or', label='Fit')
ax.loglog(kxz, Fw_exact(kxz), 'k', label='Exact')
ax.set_ylabel(r"$\mathrm{k_xF_w(k_x)}/u_{\tau}^2$")
ax.set_xlabel(r"$\mathrm{fz}/U_0$")
savefig('w_spectrum.png')
EOF
else
    exit 1
fi
