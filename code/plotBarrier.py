import numpy as np
import matplotlib.pyplot as plt

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

x,psiR,psiI,psi2 = np.genfromtxt("../data/pot_barrier_wave_k0_200.dat",unpack=True,skip_header=True)

plt.figure("Potential Well")
plt.plot(x,psi2,label="Numerical solution")
#plt.plot(x,psiR,label="REal")
#plt.plot(x,psiI,label="imag")
#plt.plot(xx,abs(psiEnd)**2,'--',label="Analytic solution")
plt.plot(xx,psi02,label="Initial wavefunction")
plt.plot(xx,V,label="Potential")
#plt.xlim(-4,4)
plt.legend()
if savePlots:
    plt.savefig("../figs/pot-barrier.pdf")


t,R,T = np.genfromtxt("../data/pot_barrier_RT_k0_150.dat",unpack=True,skip_header=True)

plt.figure("relfection transmission")
plt.plot(t,R,label="Reflection coefficient")
plt.plot(t,T,label="Transmission coefficient")
plt.legend()

k = np.linspace(150,350,100)
V0 = 2e4
E = k**2/2
kappa = np.sqrt(2*E)

Tana = 1/(1+(V0**2/(4*E(V0-E))*(np.sinh(2*a*np.sqrt(2*(V0-E))))**2))
Ranal = 0*k -Tana+1

k0,t,R,T = np.genfromtxt("../data/pot_barrier_RT_all.dat",unpack=True)

plt.figure("RT vs k0")
plt.plot(k0,R,'o',label="Reflection")
plt.plot(k0,T,'o',label="Transmission")
plt.plot(k,Tana,label="aTransmission")
plt.plot(k,Ranal,label="aReflection")
plt.axvline(200)
plt.legend()

if showPlots:
    plt.show()
