import numpy as np
import matplotlib.pyplot as plt

savePlots = 0
showPlots = 1

x0=-2;
s0 = 0.1;
t1 = 8e-3;
xx = np.linspace(-8,8,49152)
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
plt.xlim(-3,3)
plt.ylabel(r"$|\Psi|^{2}$")
plt.xlabel("x")
plt.legend()
if savePlots:
    plt.savefig("../figs/free-part-norm2.pdf")

plt.figure("imag part")
plt.plot(x,psiI,label="Numerical solution")
plt.plot(xx,psiEnd.imag,'--',label="Analytic solution")
plt.xlim(-0.5,0.35)
plt.ylabel(r"$\mathfrak{Im}(\Psi)$")
plt.xlabel("x")
plt.legend()
if savePlots:
    plt.savefig("../figs/free-part-im.pdf")

if showPlots:
    plt.show()

print(np.max(np.abs(psi2-abs(psiEnd)**2)))
