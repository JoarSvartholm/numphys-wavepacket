import numpy as np
import matplotlib.pyplot as plt

x0=-2;
s0 = 0.1;
t1 = 8e-3;
xx = np.linspace(-8,8,10000)
psi02 = np.sqrt(1/np.pi)/s0*np.exp(-(xx-x0)**2/(s0**2))

k0 = 240;
theta = 0.5*np.arctan(t1/s0)
phi = -theta - k0**2*t1/2
psiEnd =(s0**2/np.pi)**0.25 * np.exp(1j*phi)/(s0**2+1j*t1)**0.5 * np.exp(1j*k0*xx) * np.exp( -(xx-x0-k0*t1)**2/(2*s0**2+2*1j*t1) )


x,psiR,psiI,psi2 = np.genfromtxt("../data/free_part_wave",unpack=True,skip_header=True)

plt.figure("free particle")
plt.plot(x,psi2,label="Numerical solution")
plt.plot(xx,abs(psiEnd)**2,'--',label="Analytic solution")
plt.plot(xx,psi02,label="Initial wavefunction")
plt.xlim(-4,4)
plt.legend()
plt.show()
