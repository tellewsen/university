# importing useful libraries
from functions import *
import matplotlib.pyplot as plt  # plotting package
from matplotlib import rc 
rc('font', **{'family': 'serif'})  # This is for Latex writing

# Part3
# 3.1 The Planck Law
wav = np.arange(1000, 40001, 200)
b = np.zeros(wav.shape)

plt.figure(0)
plt.xlabel(r'wavelength $\lambda / \AA$', size=14)
plt.ylabel(r'Planck function', size=14)
plt.title('Radiation intensity emitted')
plt.set_cmap('RdYlBu')
for T in np.arange(5000, 8001, 500):
    b[:] = planck(T, wav[:]*1e-8)
    plt.plot(wav, b, '-')

plt.yscale('log')
plt.xscale('log')
plt.xlim(0, 40001)

# 3.2 Radiation through an isothermal layer

B = 2.
tau = np.arange(0.01, 10, 0.01)
intensity = np.zeros(tau.shape)

plt.figure(1)
for I_0 in range(4, -1, -1):
    intensity[:] = I_0*np.exp(-tau[:]) + B*(1.-np.exp(-tau[:]))
    plt.plot(tau, intensity, label=r'intensity I$_0$ = ' + str(I_0))

plt.xlabel(r'optical depth $\tau$', size=14)
plt.ylabel(r'intensity', size=14)
plt.legend(loc='best', fontsize=12)
plt.title('Emergent intensity through a layer')
plt.xscale('log')
plt.yscale('log')


# 3.3 Spectral lines from a solar reversing layer
u = np.arange(-10, 10.1, 0.1)
a = np.array([0.001, 0.01, 0.1, 1])
vau = np.zeros((a.shape[0], u.shape[0]))

plt.figure(2)
for i in range(4):
    vau[i, :] = voigt(a[i], u[:])
    plt.plot(u[:], vau[i, :], label='a = ' + np.str(a[i]))
plt.title('Voigt profile for different a values')
plt.ylim(0, 1)
plt.xlim(-10, 10)
plt.legend(fontsize=12)
plt.xlabel('u')
plt.ylabel('voigt profile', size=12)

plt.figure(3)
plt.title('Voigt profile for different a values')
for i in range(4):
    vau[i, :] = voigt(a[i], u[:])
    plt.plot(u[:], vau[i, :], label='a = ' + np.str(a[i]))
plt.yscale('log')
plt.legend(fontsize=12, loc=8)
plt.xlabel('u', size=14)
plt.ylabel('logarithmic voigt profile', size=12)
plt.show()

# 3.3 Emergent line profiles
Ts = 4200.  # Temperature of the surface of the star
Tl = 5700.  # Temperature of the reversing layer
a = 0.01  # Coulomb damping parameter
wav = 5000e-8  # Wavelength in cm
tau_0 = 1.  # reversing layer thickness at line center
u = np.arange(-10, 10, 0.1)
intensity = np.zeros(u.shape)

for i in range(len(u)):
    tau = tau_0*voigt(a, u[i])
    intensity[i] = planck(Ts, wav)*np.exp(-tau) + planck(Tl, wav)*(1.-np.exp(-tau))
plt.figure(4)
plt.plot(u, intensity)
plt.xlabel('u')
plt.ylabel('intensity')
plt.title('Intensity around center of spectral line')


plt.figure(5)
for i_wav in range(1, 4):
    wav = (i_wav ** 2 + 1.) * 1e-5
    for i in range(len(u)):
        tau = tau_0*voigt(a, u[i])
        intensity[i] = planck(Ts, wav)*np.exp(-tau) + planck(Tl, wav)*(1.-np.exp(-tau))
    plt.plot(u, intensity, label=r'$\lambda = $' + (np.str(int(wav*1e8))) + r'$\AA$')

plt.legend()
plt.xlabel('u')
plt.ylabel('intensity')
plt.title('Intensity around center of spectral line')


log_tau_0 = np.arange(-2, 2.1, 0.5)
plt.figure(6)
for element in log_tau_0:
    for i in range(len(u)):
        tau = 10.**element * (voigt(a, u[i]))
        intensity[i] = planck(Ts, wav)*np.exp(-tau) + planck(Tl, wav)*(1.-np.exp(-tau))
    plt.plot(u, intensity, label=r'$\log{(\tau_0)} = $' + np.str(element))
    # print 'I_cont = '+np.str(intensity[0])
    # print 'I_center = ' + np.str(min(intensity))

plt.legend(loc=3, fontsize=12)
plt.xlabel('u')
plt.ylabel('intensity')
plt.title(r'Intensity as a function of u for different $\tau_0$')

plt.figure(7)
for i_wav in range(1, 4):
    wav = (i_wav ** 2 + 1.) * 1e-5
    for element in log_tau_0:
        for i in range(len(u)):
            tau = 10.**element * (voigt(a, u[i]))
            intensity[i] = planck(Ts, wav)*np.exp(-tau) + planck(Tl, wav)*(1.-np.exp(-tau))
        intensity = intensity / intensity[0]
        plt.plot(u, intensity[:], linewidth=1.)
plt.xlabel('u')
plt.ylabel(r'$I_\lambda/I_0$')
plt.title('Intensity scaled to local continuum')


# 3.4 Equivalent width of spectral lines
# Checking the profile
u = np.arange(-200, 200.4, 0.4)
a = 0.1
tau0 = 1.e2
intensity = profile(a, tau0, u, Ts, wav, Tl)

plt.figure(8)
plt.plot(u, intensity)
plt.title('Schuster-Schwarzschild profile')
plt.ylabel('Intensity')
plt.xlim([-200, 200])
plt.xlabel('u')

# relative
rel_depth = (intensity[0] - intensity) / intensity[0]

plt.figure(9)
plt.plot(u, rel_depth)
plt.title('Line depth in relative units')
plt.ylabel('Line depth')
plt.xlim([-200, 200])
plt.xlabel('u')
plt.show()

eqw = sum(rel_depth) * 0.4
print eqw

# 3.5 The curve of growth
tau0 = np.logspace(-2, 4, 61)
eqw = np.zeros(tau0.size)
u = np.arange(-200, 200.4, 0.4)
for i in range(61):
    intensity = profile(a, tau0[i], u, Ts, wav, Tl)
    rel_depth = (intensity[0] - intensity) / intensity[0]
    eqw[i] = sum(rel_depth) * 0.4

plt.figure(10)
plt.plot(tau0, abs(eqw))
plt.xlabel(r'$\tau_0$', size=18)
plt.ylabel(r'equivalent width $W_{\lambda}$', size=14)
plt.xscale('log')
plt.yscale('log')
plt.title('Curve of growth')
plt.show()
