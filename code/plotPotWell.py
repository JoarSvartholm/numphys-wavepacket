import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

savePlots = 0
showPlots = 1

x0=-1;
s0 = 0.1;
t1 = 8e-3;
xx = np.linspace(-8,8,10000)
psi02 = np.sqrt(1/np.pi)/s0*np.exp(-(xx-x0)**2/(s0**2))

k0 = 240;
theta = 0.5*np.arctan(t1/s0)
phi = -theta - k0**2*t1/2
psiEnd =(s0**2/np.pi)**0.25 * np.exp(1j*phi)/(s0**2+1j*t1)**0.5 * np.exp(1j*k0*xx) * np.exp( -(xx-x0-k0*t1)**2/(2*s0**2+2*1j*t1) )

V0 = -1
a = 0.04
V = 0*(abs(xx)>=a)+ V0*(abs(xx)<a)

x,psiR,psiI,psi2 = np.genfromtxt("../data/pot_well_wave_k0_280.dat",unpack=True,skip_header=True)

m = np.max(np.abs(psi02))
plt.figure("Potential Well")
plt.plot(x,psi2/m,label="Numerical solution")
plt.plot(xx,psi02/m,label="Initial wavefunction")
plt.plot(xx,V/m,label="Potential")
plt.xlim(-3,3)
plt.ylabel("$\Psi$ [Normalized]")
plt.xlabel("x")
plt.legend()
if savePlots:
    plt.savefig("../figs/pot-well-wave280.pdf")


t,R,T = np.genfromtxt("../data/pot_well_RT_k0_280.dat",unpack=True,skip_header=True)

plt.figure("relfection transmission")
plt.plot(t,R,label="$R$")
plt.plot(t,T,label="$T$")
plt.xlabel("t [s]")
plt.xlim(0,0.008)
plt.legend()
if savePlots:
    plt.savefig("../figs/pot-well-RT280.pdf")

k = np.linspace(150,350,100)
V0 = 1e5
E = k**2/2

Tana = 1/(1+(V0**2/(4*E*(E+V0)))*(np.sin(2*a*np.sqrt(2*(E+V0))))**2)
Ranal = 0*k -Tana+1

n = np.array((13,14))
Eres = (n*np.pi/(2*a))**2/2-V0
res = np.sqrt(2*Eres)
print(res)

k0,t,R,T = np.genfromtxt("../data/pot_well_RT_all.dat",unpack=True)

Rint = interp1d(k0,R,kind="cubic")
Tint = interp1d(k0,T,kind="cubic")

plt.figure("RT vs k0")
plt.plot(k,Rint(k),'b-',label="$R_{num}$")
plt.plot(k0,R,'bo')
plt.plot(k0,T,'ro')
plt.plot(k,Tint(k),'r-',label="$T_{num}$")
plt.plot(k,Tana,label="$T_{analytic}$")
plt.plot(k,Ranal,label="$R_{analytic}$")
plt.axvline(res[0])
plt.axvline(res[1],label="$k_{resonance}$")
plt.legend()
plt.xlabel("$k_0$ [1/m]")
if savePlots:
    plt.savefig("../figs/pot-well-RTk0.pdf")

if showPlots:
    plt.show()
