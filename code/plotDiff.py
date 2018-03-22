import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

savePlots = 1
showPlots = 0

x0=-1;
s0 = 0.1;
t1 = 8e-3;
xx = np.linspace(-8,8,10000)
psi02 = np.sqrt(1/np.pi)/s0*np.exp(-(xx-x0)**2/(s0**2))

k0 = 240;
theta = 0.5*np.arctan(t1/s0)
phi = -theta - k0**2*t1/2
psiEnd =(s0**2/np.pi)**0.25 * np.exp(1j*phi)/(s0**2+1j*t1)**0.5 * np.exp(1j*k0*xx) * np.exp( -(xx-x0-k0*t1)**2/(2*s0**2+2*1j*t1) )

V0 = 1
a = 0.01
V = 0*(abs(xx)>=a)+ V0*(abs(xx)<a)

x,psiR,psiI,psi2 = np.genfromtxt("../data/diff_barrier_wave_k0_170.dat",unpack=True,skip_header=True)

m = np.max(np.abs(psi02))
plt.figure("Potential Well")
plt.plot(x,psi2/m,label="Numerical solution")
plt.plot(xx,psi02/m,label="Initial wavefunction")
plt.plot(xx,V/m,label="Potential")
plt.xlim(-3,3)
plt.xlabel("x")
plt.ylabel("$\Psi$ [Normalized]")
plt.legend()
if savePlots:
    plt.savefig("../figs/diff-barrier-wave170.pdf")


t,R,T = np.genfromtxt("../data/diff_barrier_RT_k0_170.dat",unpack=True,skip_header=True)

plt.figure("relfection transmission")
plt.plot(t,R,label="$R$")
plt.plot(t,T,label="$T$")
plt.xlabel("t [s]")
plt.xlim(0,0.01)
plt.legend()
if savePlots:
    plt.savefig("../figs/diff-barrier-RT170.pdf")

kk= np.linspace(150,350,100)
k = np.linspace(150,199.99,100)
V0 = 2e4
E = k**2/2
kappa = np.sqrt(2*(V0-E))

Tana = 1/(1+(V0**2/(4*E*(V0-E))*(np.sinh(2*a*np.sqrt(2*(V0-E))))**2))
#Tana = 1/(1+(kappa**2+k**2)**2/(4*k**2*kappa**2)*(np.sinh(2*kappa*a)))
Ranal = 0*k -Tana+1

k2 = np.linspace(200.01,350,100)
E = k2**2/2
Tana2 = 1/(1+(V0**2/(4*E*(E-V0)))*(np.sin(2*a*np.sqrt(2*(E-V0))))**2)
Ranal2 = 0*k -Tana2+1

k0,t,R,T = np.genfromtxt("../data/diff_RT_all.dat",unpack=True)

Rint = interp1d(k0,R,kind="cubic")
Tint = interp1d(k0,T,kind="cubic")

plt.figure("RT vs k0")
plt.plot(kk,Rint(kk),'b-',label="$R_{num}$")
plt.plot(k0,R,'bo')
plt.plot(k0,T,'ro')
plt.plot(kk,Tint(kk),'r-',label="$T_{num}$")
plt.plot(k,Tana,color=(0.3,0.6,0.8),label="$T_{analytic}$")
plt.plot(k2,Tana2,color=(0.3,0.6,0.8))
plt.plot(k,Ranal,color=(0.1,0.7,0),label="$R_{analytic}$")
plt.plot(k2,Ranal2,color=(0.1,0.7,0))
plt.axvline(200,label="$k_V$")
plt.xlabel("$k_0$ [1/m]")
plt.legend()
if savePlots:
    plt.savefig("../figs/diff-barrier-RTk0.pdf")

if showPlots:
    plt.show()
