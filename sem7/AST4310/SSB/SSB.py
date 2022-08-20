# importing useful libraries
import matplotlib.pyplot as plt  # plotting package

from functions import *
from matplotlib import rc
rc('font', **{'family': 'serif', 'size': 14})  # This is for Latex writing


# Definition of constants and variables
keV = 8.61734E-5 	# Boltzmann constant in eV/K
kerg = 1.380658E-16	 # Boltzmann constant in erg K
kJoule = kerg*1e-7	 # Boltzmann constant in joule/K
h = 6.62607E-27	 # Planck constant in erg s
c = 2.99792E10  # Speed of light in cm/s
m_e = 9.109390E-28  # Electron mass in grams
m_H = 1.67352E-24   # Hydrogen mass in grams
m_He = 3.98*m_H 	# Helium mass in grams
m_p = 1.67262E-24  # Proton mass in grams
m_u = 1.660538e-24  # Atomic mass unit in grams
m_Na = 22.99*1.6605*1e-24  # Sodium mass in grams
R_sun = 696300  # Radius sun in km
D_sun = 1.496e8 	# Distance sun-earth in km
erg2eV = 6.242e11  # Convert erg to ev

# Part 1
# Read in FALC data and arrange in arrays
falc = np.loadtxt('falc.dat', unpack=True)
height = np.array(falc[0, :])
tau5 = np.array(falc[1, :])
colm = np.array(falc[2, :])
temp = np.array(falc[3, :])
vturb = np.array(falc[4, :])
nhyd = np.array(falc[5, :])
nprot = np.array(falc[6, :])
nel = np.array(falc[7, :])
ptot = np.array(falc[8, :])
pgasptot = np.array(falc[9, :])
dens = np.array(falc[10, :])

# 1.1
# Plot falc model
plt.figure(0)
plt.plot(height, temp)
plt.ylim([3000, 10000])
plt.title(r'FALC Model')
plt.xlabel(r'Height [km]', size=14)
plt.ylabel(r'Temperature [K]', size=14)

# 1.2
# Plot total pressure vs column mass linearly
plt.figure(1)
plt.plot(colm, ptot*1e-4)
plt.title(r'$P_{total}$ vs column mass m')
plt.ylabel(r'$P_{total}$ [$10^4$dyn cm$^{-2}$]', size=14)
plt.xlabel(r'Column mass m [$g cm^{-2}$]', size=14)

# Logarithmic
plt.figure(2)
plt.plot(colm, ptot)
plt.title(r'$P_{total}$ vs column mass m logarithmically')
plt.ylabel(r'log$P_{total}$ [dyn cm$^{-2}$]', size=14)
plt.xlabel(r'Column mass m [g cm$^{-2}$]', size=14)
plt.yscale('log')
# plt.xscale('log')


# Find average c value
g_S = np.average(ptot/colm)
print 'g_S = ', g_S

# Plot ratio hydrogen mass density/total mass density vs height
rho_H = nhyd * m_H	 # Hydrogen mass density
nhel = nhyd*.1  # Helium number density
rho_He = nhel*m_He  # Helum mass density

plt.figure(3)
plt.plot(height, rho_H/dens)
plt.plot(height, rho_He/dens)
plt.plot(height, (rho_H+rho_He)/dens)
plt.legend([r'$\rho_{H}/\rho_{total}$',
            r'$\rho_{He}/\rho_{total}$',
            r'$(\rho_{H}+\rho_{He})/\rho_{total}$'], loc='best')
plt.title(r'Density fractions of total as a function of height')
plt.ylim([0, 1.1])
plt.ylabel(r'$\rho/\rho_{total}$', size=14)
plt.xlabel(r'Height [km]', size=14)

# Metal fraction
metalfraction = 1. - np.average((rho_H+rho_He)/dens)
print 'Metal fraction = ', metalfraction

# Plot column mass vs height
plt.figure(4)
plt.plot(height, colm)
# plt.legend([])
plt.title(r'Column mass vs height')
# plt.ylim([0,1.1])
plt.ylabel(r'm [g cm$^{-2}$]', size=14)
plt.xlabel(r'Height [km]', size=14)


# Logarithmic
plt.figure(5)
plt.plot(height, colm)
# plt.legend([])
plt.title(r'Column mass vs height logarithmic')
# plt.ylim([0,1.1])
plt.yscale('log')
plt.ylabel(r'log(m) [g cm$^{-2}$]', size=14)
plt.xlabel(r'Height [km]', size=14)


# Scale height
rhodive = dens[-1]/np.exp(1)
rhoe = np.zeros(np.size(height))
rhoe.fill(rhodive)

# Plot gas density vs height

plt.figure(6)
plt.plot(height, dens)
plt.plot(height, rhoe)
# plt.legend([])
plt.title(r'Gas density vs height')
# plt.ylim([0,1.1])
# plt.yscale('log')
plt.ylabel(r'Gas density $\rho$ [g cm$^{-3}$]', size=14)
plt.xlabel(r'Height [km]', size=14)
plt.legend(['Gas density', r'$\rho/e$'])


# Compute gas pressure and plot against height
pgas = pgasptot*ptot  # Gas pressure
idealgas = (nhyd + nel)*kerg*temp

# Over plot product (nhyd + nel)kT
plt.figure(7)
plt.plot(height, pgas*1e-4)  # from falc
plt.plot(height, idealgas*1e-4)  # from idealgas law
plt.legend(['Falc', 'Ideal gas'])
plt.title(r'Gas pressure vs height')
# plt.ylim([0,1.1])
# plt.yscale('log')
plt.ylabel(r'Gas pressure $P_{gas}$ [$10^4$ dyn cm$^{-2}$]', size=14)
plt.xlabel(r'Height [km]', size=14)


# Plot ratio between curves

plt.figure(8)
plt.plot(height, pgas/idealgas)  # ratio curves
plt.title(r'Ratio gas pressure from falc and ideal gas law')
plt.ylabel(r'Gas pressure $P_{gas}$ [dyn cm$^{-2}$]', size=14)
plt.xlabel(r'Height [km]', size=14)

# Add n_helium and plot again
idealgas_hel = (nhyd + nel + nhel)*kerg*temp

# Overplot product (nhyd + nel + nhel)kT

plt.figure(9)
plt.plot(height, pgas*1e-4)  # from falc
plt.plot(height, idealgas_hel*1e-4)  # from ideal gas with helium
plt.legend(['Falc', 'Ideal gas with helium'])
plt.title(r'Gas pressure vs height')
# plt.ylim([0,1.1])
# plt.yscale('log')
plt.ylabel(r'$P_{gas}$ [$10^4$ dyn cm$^{-2}$]', size=14)
plt.xlabel(r'Height [km]', size=14)


# Plot ratio between curves

plt.figure(10)
plt.plot(height, pgas/idealgas_hel)  # ratio curves
plt.title(r'Ratio gas pressure from falc and ideal gas law')
plt.ylabel(r'Ratio $P_{FALC}/P_{IG}$', size=14)
plt.xlabel(r'Height [km]', size=14)
plt.ylim([0, 1.1])


# Plot total hydrogen density vs height and overplot electron density,
# proton density, and density of electrons not from hydrogen ionization.

n_eb = (nhyd - nprot)

plt.figure(11)
plt.plot(height, nhyd)  # hydrogen vs h
plt.plot(height, nel)  # electron vs h
plt.plot(height, nprot)  # proton vs h
plt.plot(height, n_eb)  # bound electron vs h
plt.legend([r'$n_H$', r'$n_e$', r'$n_p$', r'$n_{be}$'])
plt.title(r'Number Densities vs height')
plt.ylabel(r'$n$ [cm$^{-3}$]', size=14)
plt.xlabel(r'Height [km]', size=14)
# plt.yscale('log')


# Plot ionization fraction of hydrogen logarithmically against height
hyd_ion = nprot/nhyd

plt.figure(12)
plt.plot(height, hyd_ion)
plt.title(r'Ionization fraction of hydrogen vs height')
plt.ylabel(r'log($n_p/n_H$)', size=14)
plt.xlabel(r'Height [km]', size=14)
plt.yscale('log')
plt.xlim([-100, 2220])


# Photon density
N_phot = 20*temp**3

print 'N_phot for deep photosphere = ', N_phot[np.where(height == np.min(height))][0]
print 'N_H for deep photosphere    = ', nhyd[np.where(height == np.min(height))][0]
Teff = 5770.
N_phot_high = 20.*Teff**3/(2.*np.pi)
print 'N_phot for highest point sun = %e' % N_phot_high
print 'N_H for highest point sun    = %e' % nhyd[np.where(height == np.max(height))][0]

# Write code to read earth.dat
earth = np.loadtxt('earth.dat', unpack=True)
eh = np.array(earth[0, :])
elogP = np.array(earth[1, :])
etemp = np.array(earth[2, :])
elogdens = np.array(earth[3, :])
elogN = np.array(earth[4, :])

# Plot temp,pressure,density,gas density

plt.figure(13)
plt.plot(eh, elogP)  # ratio curves
plt.title(r'Air pressure vs height logarithmic')
plt.ylabel(r'log$P$ [dyn cm$^{-2}$]', size=14)
plt.xlabel(r'Height [km]', size=14)


plt.figure(14)
plt.plot(eh, etemp)  # ratio curves
plt.title(r'Temperature vs height')
plt.ylabel(r'Temperature [K]', size=14)
plt.xlabel(r'Height [km]', size=14)


plt.figure(15)
plt.plot(eh, elogdens)  # ratio curves
plt.title(r'Gas density vs height logarithmic')
plt.ylabel(r'log $\rho$ [g cm$^{-3}$]', size=14)
plt.xlabel(r'Height [km]', size=14)


plt.figure(16)
plt.plot(eh, elogN)  # ratio curves
plt.title(r'Particle density vs height')
plt.ylabel(r'log $N$ [cm$^{-3}$]', size=14)
plt.xlabel(r'Height [km]', size=14)

# Plot pressure and density stratifications together in normalized units in one graph
eP = 10.**elogP	 # Convert from log to actual
edens = 10.**elogdens  # Convert from log to actual
eN = 10.**elogN
edense = np.zeros(np.size(eh))
edense.fill(edens[np.where(eh == 0)][0]/np.exp(1))

plt.figure(17)
plt.plot(eh, edens)
plt.plot(eh, edense)

plt.figure(18)
plt.plot(eh, edens/np.max(edens))
plt.plot(eh, eP/np.max(eP))
plt.plot(eh, eN/np.max(eN))
plt.yscale('log')  # plot on logscale
plt.legend([r'$\rho/\rho_{max}$', '$P/P_{max}$', '$N/N_{max}$'])
plt.title(r'Density and pressure stratification')
plt.ylabel(r'Values')
plt.xlabel(r'Height [km]')

# Plot mean molecular weight vs height
my_E = edens/(eN*m_H)  # mean molecular weight

plt.figure(19)
plt.plot(eh, my_E)
plt.title(r'Mean molecular weight vs height')
plt.ylabel(r'$\mu_E$[]')
plt.xlabel(r'Height [km]')
plt.show()

# Estimate density scale height of lower terrestrial atmosphere
g_E = 980.665  # cm s^-2
eH_p = kJoule*etemp/((my_E*m_u*1e-3)*(g_E*1e-2))*1e-3
print 'Scale height H_p low earth atmosphere =  ', eH_p[0], ' km'

# Compare terrestrial parameter values to solar
print '------------------'
print 'All values calculated at h = 0'

print 'Matter density sun   = ', dens[np.where(height == 0)][0]
print 'Matter density earth = ', edens[np.where(eh == 0)][0]

N = (nhyd + nhel + nel)
print 'Particle density sun   = ', N[np.where(height == 0)][0]
print 'Particle density earth = ', eN[np.where(eh == 0)][0]

eptot = 10**elogP

print 'Pressure sun   = ', ptot[np.where(height == 0)][0]
print 'Pressure earth = ', eptot[np.where(eh == 0)][0]

print 'Temperature sun   = ', temp[np.where(height == 0)][0]
print 'Temperature earth = ', etemp[np.where(eh == 0)][0]

print 'Ratio particle densities = ', eN[np.where(eh == 0)][0]/N[np.where(height == 0)][0]
print '------------------'

# Estimate atmospheric column mass earth
g_E = 980.665  # cm s^-2
ecolm = eptot[np.where(eh == 0)][0]/g_E
print 'Column mass earth surface = ', ecolm
print 'Column mass sun surface = ', colm[np.where(height == 0)][0]

# Compare N_phot to particle density in the air on earth
N_phot_sun = np.pi*R_sun**2/D_sun**2*N_phot_high
N_phot_earth = 20*etemp[np.where(eh == 0)][0]**3
print 'Photons from sun: %e' % N_phot_sun
print 'Photons from earth:%e' % N_phot_earth
print 'Particle density surface of earth: %e' % eN[np.where(eh == 0)][0]
print 'Ratio photons from sun and particles earth atmosphere = ', N_phot_sun/eN[np.where(eh == 0)][0]
print 'Ratio photons from sun and photons from earth atmosphere = ', N_phot_sun/N_phot_earth


# Part 2 Observed solar continua

# Write code to read table 5 (solspect.dat)
solspect = np.loadtxt('solspect.dat', unpack=True)
wavelength = solspect[0, :]

F_lambda = solspect[1, :]
F_lambda_c = solspect[2, :]
I_lambda = solspect[3, :]
I_lambda_c = solspect[4, :]

# Plot the four distributions in one figure over the range wavelength = 0-2 micro meter.

plt.figure(20)
plt.plot(wavelength, F_lambda)
plt.plot(wavelength, F_lambda_c)
plt.plot(wavelength, I_lambda)
plt.plot(wavelength, I_lambda_c)
plt.xlim([0, 2])
plt.title('The four spectral distributions')
plt.legend([r'$F_\lambda$', r'$F_\lambda^c$', r'$I_\lambda$', r'$I_\lambda^c$'])
plt.xlabel(r'Wavelength $\lambda$ [$\mu m$]')
plt.ylabel(r'Distributions [$10^{10}$erg cm$^{-2}$s$^{-1}$ster$^{-1}\mu$m$^{-1}$]')


# Check that continuum intensity reaches 4.6*10^10 at max
print 'I_lambda_c = %e' % np.max(I_lambda_c)  # check okay

# Convert distributions into values per frequency
conversionfactor = 1./c*(wavelength*1e-4)**2*1e4
F_nu = F_lambda*conversionfactor
F_nu_c = F_lambda_c*conversionfactor
I_nu = I_lambda*conversionfactor
I_nu_c = I_lambda_c*conversionfactor

frequency = c/(wavelength*1e-4)

plt.figure(21)
plt.plot(wavelength, F_nu*1e15)
plt.plot(wavelength, F_nu_c*1e15)
plt.plot(wavelength, I_nu*1e15)
plt.plot(wavelength, I_nu_c*1e15)
plt.title('The four spectral distributions')
plt.legend([r'$F_\nu$', r'$F_\nu^c$', r'$I_\nu$', r'$I_\nu^c$'])
plt.xlabel(r'Wavelength $\lambda$ [$\mu m$]')
plt.ylabel(r'Distributions [$10^{-5}$erg cm$^{-2}$s$^{-1}$ster$^{-1}$Hz$^{-1}$]')
plt.xlim([0, 2])


print 'Peak I_nu_c = %e' % np.max(I_nu_c)

planckfit1 = planck(6200, wavelength)
planckfit2 = planck(6300, wavelength)
planckfit3 = planck(6400, wavelength)

plt.figure(22)
plt.plot(wavelength, I_lambda_c)
plt.plot(wavelength, planckfit1)
plt.plot(wavelength, planckfit2)
plt.plot(wavelength, planckfit3)
plt.xlim([0, 2])
plt.title(r'Observed intensity of continuum vs Planck fit')
plt.ylabel(r'Intensity [$10^{10}$ erg cm$^{-2}$s$^{-1}$ster$^{-1}\mu$m$^{-1}$]')
plt.xlabel(r'Wavelength $\lambda$ [$\mu m$]')
plt.legend([r'$I_\lambda^c$', r'$B_\lambda(6200)$', r'$B_\lambda(6300)$', r'$B_\lambda(6400)$'])


# planckfit for frequency
planckfit4 = planck_nu(6200, wavelength)*1e15  # 1e15 changes units on y axis.
planckfit5 = planck_nu(6300, wavelength)*1e15  # note the change from 10^10
planckfit6 = planck_nu(6400, wavelength)*1e15  # to 10^-5.

plt.figure(23)
plt.plot(wavelength, I_nu_c*1e15)
plt.plot(wavelength, planckfit4)
plt.plot(wavelength, planckfit5)
plt.plot(wavelength, planckfit6)
plt.xlim([0, 2])
plt.title(r'Observed intensity of continuum vs Planck fit')
plt.ylabel(r'Intensity [$10^{-5}$ erg cm$^{-2}$s$^{-1}$ster$^{-1}$Hz$^{-1}$]')
plt.xlabel(r'Wavelength $\lambda$ [$\mu m$]')
plt.legend([r'$I_\nu^c$', r'$B_\nu(6200)$', r'$B_\nu(6300)$', r'$B_\nu(6400)$'])


T_b = planck_invert(I_lambda_c, wavelength)  # Brightness temperature solar continuum

plt.figure(24)
plt.plot(wavelength, T_b)
plt.title(r'Brightness temperature vs wavelength for the solar continuum')
plt.ylabel(r'Brightness temperature $T_b$[K] ')
plt.xlabel(r'Wavelength [$\mu$m]')


# 2.2 Continuous extinction
wav = np.linspace(3000, 19000, 1000)
FALCexthmin = exthmin(wav, temp[np.where(height == 0)][0], nel[np.where(height == 0)][0])
print temp[np.where(height == 0)][0]
print nel[np.where(height == 0)][0]
print r'Maximum of exhtmin = %e at wavelength = %s angstrom' % \
      (np.max(FALCexthmin), wav[np.where(FALCexthmin == np.max(FALCexthmin))][0])

plt.figure(25)
plt.title(r'H-minus bound-free + free-free extinction per H atom')
plt.plot(wav*1e-4, FALCexthmin*1e26/(nel[np.where(height == 0)][0] *
                                     kerg*temp[np.where(height == 0)][0]))  # Note the 10^-24 units
plt.ylabel(r'$\kappa /P_e [10^{-26}$cm$^2$/dyn cm$^{-2}$]', size=14)
plt.xlabel(r'Wavelength [$\mu$m]')
plt.xlim([0, 2])
plt.ylim([0, 5])
plt.xscale('log')
plt.yscale('log')


# Plot variation of H^- extinction with height for wav = 0.5 \mu m
exthmin05 = exthmin(5000, temp, nel)*(nhyd - nprot)
thomson = 6.648*1e-25  # cm**2
extincelec = thomson*nel
exttotal = exthmin05+extincelec

plt.figure(26)
plt.plot(height, exthmin05)
plt.plot(height, extincelec)
plt.plot(height, exttotal)
plt.legend([r'Neutral H', r'Free electrons', r'Total'])
plt.title(r'H-minus bf+ff extinction at $\lambda = 0.5 \mu m$')
plt.ylabel(r'log $(\kappa)$ [cm$^{-1}$]', size=14)
plt.xlabel(r'Height [km]')
plt.yscale('log')

# part 2.3
# integrate the extinction
tau = np.zeros(np.size(tau5))  # initializing tau array
for i in range(1, len(tau)):  # index zero is not accounted for, so tau[0] = 0 because we have air
    tau[i] = tau[i-1] + 0.5*(exttotal[i]+exttotal[i-1])*(height[i-1]-height[i])*1e5

plt.figure(27)
plt.plot(height, tau5, '--')
plt.plot(height, tau)
plt.yscale('log')
plt.legend(['FALC', 'Integrated'])
plt.title(r'$\tau_{500}$ from FALC compared to numerical integration')
plt.ylabel(r'Optical depth $\tau_{500}$')
plt.xlabel(r'Height [km]')

# 2.5 Emergent intensity
ext = np.zeros(np.size(tau5)) 
tau = np.zeros(np.size(tau5))
integrand = np.zeros(np.size(tau5)) 
contfunc = np.zeros(np.size(tau5))  # repeat for every parameter if e
intt = 0.0 
hint = 0.0
wl = 5.  # wavelength in micrometers
for i in range(1, len(tau)):  # the index zero is not accounted for
    ext[i] = exthmin(wl*1e4, temp[i], nel[i])*(nhyd[i]-nprot[i])+0.664e-24*nel[i]
    tau[i] = tau[i-1] + 0.5*(ext[i]+ext[i-1])*(height[i-1]-height[i])*1e5
    integrand[i] = planck(temp[i], wl)*np.exp(-tau[i])
    intt += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
    hint += height[i]*0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
    contfunc[i] = integrand[i]*ext[i]

hmean = hint / intt

# Compare the computed intensity at 500 nm with the observed one
print 'I_lambda_cont observed = ', I_lambda_c[np.where(wavelength == wl)][0]
print 'I_lambda_cont computed = ', intt
# print (1- I_lambda_c[np.where(wavelength == wl)][0]/intt)*100. #  deviation from observed
print 'hmean = ', hmean  # mean height of formation

# plot peak normalized contrib func vs height
print 'peak location = ', height[np.where(contfunc == np.max(contfunc))][0]

plt.figure(28)
plt.plot(height, contfunc/np.max(contfunc))
plt.title(r'Contribution function vs height for $\lambda = 50000\AA$')
plt.ylabel(r'$I_\lambda/I_\lambda^{peak}$', size=14)
plt.xlabel(r'Height[km]')
plt.xlim([np.min(height), 1500])


print 'tau= ', find_nearest(tau, 1)[0], 'at h= ', height[find_nearest(tau, 1)[1]]

# Print height where T is closest to T_b
print 'For wav =  500nm, T is closest to T_b at h= ', height[find_nearest(temp, T_b[np.where(wavelength == 0.5)][0])[1]]
print 'At that point T_B = ', T_b[np.where(wavelength == 0.5)][0], 'and T = ', \
    temp[find_nearest(temp, T_b[np.where(wavelength == 0.5)])[1]]

print 'For wav = 1000nm, T is closest to T_b at h= ', height[find_nearest(temp, T_b[np.where(wavelength == 1.0)][0])[1]]
print 'At that point T_B = ', T_b[np.where(wavelength == 1.0)][0], 'and T = ', \
    temp[find_nearest(temp, T_b[np.where(wavelength == 1.0)])[1]]

print 'For wav = 1600nm, T is closest to T_b at h= ', height[find_nearest(temp, T_b[np.where(wavelength == 1.6)][0])[1]]
print 'At that point T_B = ', T_b[np.where(wavelength == 1.6)][0], 'and T = ', \
    temp[find_nearest(temp, T_b[np.where(wavelength == 1.6)])[1]]

print 'For wav = 5000nm, T is closest to T_b at h= ', height[find_nearest(temp, T_b[np.where(wavelength == 5.0)][0])[1]]
print 'At that point T_B = ', T_b[np.where(wavelength == 5.0)][0], 'and T = ', \
    temp[find_nearest(temp, T_b[np.where(wavelength == 5.0)])[1]]


# 2.5 Disk center intensity

ext = np.zeros(np.size(tau5)) 
tau = np.zeros(np.size(tau5))
integrand = np.zeros(np.size(tau5)) 
contfunc = np.zeros(np.size(tau5))

intt = np.zeros(np.size(wavelength))
hint = np.zeros(np.size(wavelength))
for j in range(np.size(wavelength)):
    for i in range(1, len(tau)):  # the index zero is not accounted for
        ext[i] = exthmin(wavelength[j]*1e4, temp[i], nel[i])*(nhyd[i]-nprot[i])+0.664e-24*nel[i]
        tau[i] = tau[i-1] + 0.5*(ext[i]+ext[i-1])*(height[i-1]-height[i])*1e5
        integrand[i] = planck(temp[i], wavelength[j])*np.exp(-tau[i])
        intt[j] += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
        hint[j] += height[i]*0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
        contfunc[i] = integrand[i]*ext[i]


plt.figure(29)
plt.plot(wavelength, intt)
plt.plot(wavelength, I_lambda_c)
plt.title(r'Disk-center intensity', size=14)
plt.xlabel(r'Wavelength [$\mu$m]', size=14)
plt.ylabel(r'Intensity $I_\lambda$ [$10^{10}$ erg cm$^{-2}$s$^{-1}$ster$^{-1}\mu$m$^{-1}$]', size=14)
plt.legend([r"Computed", "Observed"])


# 2.5 Disk center intensity continued

ext = np.zeros(np.size(tau5)) 
tau = np.zeros(np.size(tau5))
integrand = np.zeros(np.size(tau5)) 
mu = np.linspace(0.1, 1, 50)
intt = np.zeros((np.size(mu), np.size(wavelength)))
hint = np.zeros((np.size(mu), np.size(wavelength)))
for k in range(0, np.size(mu)):
    for j in range(np.size(wavelength)):
        for i in range(1, len(tau)):  # the index zero is not accounted for
            ext[i] = exthmin(wavelength[j]*1e4, temp[i], nel[i])*(nhyd[i]-nprot[i])+0.664e-24*nel[i]
            tau[i] = tau[i-1] + 0.5*(ext[i]+ext[i-1])*(height[i-1]-height[i])*1e5
            integrand[i] = planck(temp[i], wavelength[j])*np.exp(-tau[i]/mu[k])/mu[k]
            intt[k, j] += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
            hint[k, j] += height[i]*0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])


# plt.figure(29)
# for k in range(np.size(mu)):
#     plt.plot(wavelength, intt[k], label=r'$\mu=%s$' % mu[k])
# plt.title(r'Disk-center intensity', size=14)
# plt.xlabel(r'Wavelength [$\mu$m]', size=14)
# plt.ylabel(r'Intensity $I_\lambda$ [$10^{10}$ erg cm$^{-2}$s$^{-1}$ster$^{-1}\mu$m$^{-1}$]', size=14)
# # plt.ylim([-.1, 6.])
# plt.legend(loc='best')
# # plt.yscale('log')


# 2.6 Limb darkening

ratio1 = np.zeros(np.size(mu))
ratio2 = np.zeros(np.size(mu))
ratio3 = np.zeros(np.size(mu))
ratio4 = np.zeros(np.size(mu))

for i in range(np.size(mu)):
    ratio1[i] = intt[i, np.where(wavelength == 0.5)]/intt[-1, np.where(wavelength == 0.5)]
    ratio2[i] = intt[i, np.where(wavelength == 1.0)]/intt[-1, np.where(wavelength == 1.0)]
    ratio3[i] = intt[i, np.where(wavelength == 1.6)]/intt[-1, np.where(wavelength == 1.6)]
    ratio4[i] = intt[i, np.where(wavelength == 5.0)]/intt[-1, np.where(wavelength == 5.0)]


plt.figure(30)
plt.plot(mu, ratio1, label=r'$\lambda = 500$nm')
plt.plot(mu, ratio2, label=r'$\lambda = 1000$nm')
plt.plot(mu, ratio3, label=r'$\lambda = 1600$nm')
plt.plot(mu, ratio4, label=r'$\lambda = 5000$nm')
plt.title(r'Ratio $I_\lambda(0,\mu)/I_\lambda(0,1)$ for select wavelengths')
plt.ylabel(r'$I_\lambda(0,\mu)/I_\lambda(0,1)$')
plt.xlabel(r'Angle $\mu$')
plt.xlim([np.max(mu), np.min(mu)])
plt.legend(loc='best')


rRsun = np.sin(np.arccos(mu))


plt.figure(31)
plt.plot(rRsun, ratio1, label=r'$\lambda = 500$nm')
plt.plot(rRsun, ratio2, label=r'$\lambda = 1000$nm')
plt.plot(rRsun, ratio3, label=r'$\lambda = 1600$nm')
plt.plot(rRsun, ratio4, label=r'$\lambda = 5000$nm')
plt.title(r'Ratio $I_\lambda(0,\mu)/I_\lambda(0,1)$ for select wavelengths')
plt.ylabel(r'$I_\lambda(0,\mu)/I_\lambda(0,1)$')
plt.xlabel(r'$r/R_\odot$')
# plt.xlim([np.max(rRsun),np.min(rRsun)])
plt.legend(loc='best')
plt.show()


# 2.7 Flux integration

# Gauss-Legendre quadrature for the integral
xgauss = [-0.7745966692, 0.0000000000, 0.7745966692]
wgauss = [0.5555555555, 0.8888888888, 0.5555555555]
fluxspec = np.zeros(np.size(wavelength))
intmu = np.zeros((3, np.size(wavelength)))
for imu in range(3):
    mu = 0.5+xgauss[imu]/2.
    wg = wgauss[imu]/2.
    for iw in range(np.size(wavelength)):
        wl = wavelength[iw]
        intt = intensity_mu(tau5, wl, mu, temp, nel, nhyd, nprot, height)[0]
        intmu[imu, iw] = intt
        fluxspec[iw] = fluxspec[iw]+wg*intmu[imu, iw]*mu
fluxspec *= 2.

plt.figure(32)
plt.plot(wavelength, fluxspec, label=r'Calculated')
plt.plot(wavelength, F_lambda_c, label=r'Observed')
plt.title(r'Solar flux')
plt.ylabel(r'$F_\lambda$ [$10^{10}$erg cm$^{-2}$s$^{-1}$ster$^{-1}\mu$m$^{-1}$]]', size=14)
plt.xlabel(r'wavelength [$\mu$m]', size=14)
plt.legend(loc='best')
plt.show()


# Part3

# Read in solar spectrum data
Nadata = np.loadtxt('int_nad.dat', unpack=True)
wavenumber = Nadata[0]
E_spec = Nadata[1]
S_spec = Nadata[2]
S_spec_corr = Nadata[3]


plt.figure(33)
plt.plot(wavenumber, S_spec)
plt.title('Observed Solar spectrum')
# plt.ylabel('')
plt.xlabel(r'Wavenumber[cm$^{-1}$]')
plt.xlim([np.min(wavenumber), np.max(wavenumber)])
# plt.show()

wavelength_vac = 1./wavenumber*1e8

plt.figure(34)
plt.plot(wavelength_vac, S_spec)
plt.title('Observed Solar spectrum')
# plt.ylabel('')
plt.xlabel(r'Wavelength[$\AA$]')
plt.xlim([np.min(wavelength_vac), np.max(wavelength_vac)])
plt.show()

# Split the spectrum into two regimes to find local minimums
S_spec1 = S_spec[:3600]
S_spec2 = S_spec[3600:]
# print wavelength of respective minima
print 'First line has minimum at : ', wavelength_vac[3600+np.argmin(S_spec2)], 'angstrom'
print 'Second line has minimum at : ', wavelength_vac[np.argmin(S_spec1)], 'angstrom'


wavelength_air = vac_to_air(wavelength_vac)
NaI_D2_minimum = vac_to_air(wavelength_vac[3600+np.argmin(S_spec2)])
NaI_D1_minimum = vac_to_air(wavelength_vac[np.argmin(S_spec1)])
print 'First line after correction has minimum at : ', NaI_D1_minimum, 'angstrom'
print 'Second line after corection has minimum at : ', NaI_D2_minimum, 'angstrom'

plt.figure(35)
plt.plot(wavelength_air, S_spec)
plt.title('Observed Solar spectrum in air')
# plt.ylabel('')
plt.xlabel(r'Wavelength[$\AA$]')
plt.xlim([np.min(wavelength_air), np.max(wavelength_air)])


# Part3.6
b_l = 1.
b_u = 1.
A_Na = 1.8*1e-6
f_lu = [0.318, 0.631]
E_ionization = np.array([5.139, 47.29, 71.64])
boltz = np.zeros((3, len(temp)))

for i in range(len(temp)):
    boltz[0, i] = boltzmann_na(temp[i], 0, 0)
    boltz[1, i] = boltzmann_na(temp[i], 0, 1)
    boltz[2, i] = boltzmann_na(temp[i], 0, 2)

plt.figure(36)
plt.plot(height, boltz[0], label=r's=1 groundstate')
plt.plot(height, boltz[1], '--', label=r's=2 upper level D1')
plt.plot(height, boltz[2], '-.', label=r's=3 upper level D2')
plt.title(r'Boltzmann distribution of Na I')
plt.xlabel(r'Height [km]')
plt.ylabel(r'Population fraction $n_{1,s}/N_1$')
plt.xlim([np.min(height), 2000])
plt.legend(loc='best')
plt.show()


# Saha
saha = np.zeros((2, len(temp)))
for i in range(len(temp)):
    saha[0, i] = saha_na(temp[i], nel[i], 1, E_ionization)
    saha[1, i] = saha_na(temp[i], nel[i], 2, E_ionization)

# plt.figure(37)
# plt.plot(height, saha[0], '-', label='Na I')
# plt.plot(height, saha[1], '--', label='Na II')
# plt.xlim([np.min(height), 2000])
# plt.ylim([1e-4, 10])
# plt.yscale('log')
# plt.title(r'Saha distribution of Na ')
# plt.xlabel(r'Height [km]')
# plt.ylabel(r'Ionization state fraction $N_r/N_{total}$')
# plt.legend(loc='best')
# plt.show()


# Saha Boltzmann

sahaboltz = np.zeros((3, len(temp)))
for i in range(len(temp)):
    sahaboltz[0, i] = sahabolt_na(temp[i], nel[i], 0, 0, E_ionization)
    sahaboltz[1, i] = sahabolt_na(temp[i], nel[i], 0, 1, E_ionization)
    sahaboltz[2, i] = sahabolt_na(temp[i], nel[i], 0, 2, E_ionization)
# plt.figure(38)
# plt.plot(height, sahaboltz[0])
# plt.plot(height, sahaboltz[1])
# plt.plot(height, sahaboltz[2])
# plt.show()


# Dopplerwidth testing

doppler_term1 = np.sqrt(2.*kerg*temp/m_Na)*1e-5
doppler_term2 = np.sqrt(vturb*vturb)

# plt.figure(39)
# plt.plot(height, doppler_term1, label=r'$\sqrt{2kT/m_{Na}}$')
# plt.plot(height, doppler_term2, '--', label=r'$v_{turb}$')
# plt.xlim(np.min(height), 2000)
# plt.ylim(0, )
# plt.ylabel('[km/s]')
# plt.xlabel('Height [km]')
# plt.title('Na Doppler broadening and micro turbulence')
# plt.legend(loc='best')
# plt.show()

doppler = np.zeros((2, len(temp)))
wav = np.array([NaI_D1_minimum, NaI_D2_minimum])*1e-8  # wavelengths of line in cm
doppler[0] = doppler_width(wav[0], temp, vturb, m_Na)  # values for NaID1
doppler[1] = doppler_width(wav[1], temp, vturb, m_Na)  # for NaID2

# plt.figure(40)
# plt.plot(height, doppler[0]*1e8, label='Na I D1')
# plt.plot(height, doppler[1]*1e8, 'r--', label='Na I D2')
# plt.xlim(np.min(height), 2000)
# plt.ylabel(r'Doppler width [$\AA$]')
# plt.xlabel(r'Height [km]')
# plt.title(r'Doppler width')
# plt.legend()
# plt.show()


# Van der Waal broadening Na I D1 line
pgas = pgasptot*ptot
vdw_broadening = np.zeros(len(temp))
for i in range(len(temp)):
    vdw_broadening[i] = gammavdw_nad(temp[i], pgas[i], 2, wav, E_ionization)

# plt.figure(41)
# plt.plot(height, vdw_broadening)
# plt.yscale('log')
# plt.xlim([np.min(height), 2000])
# plt.title('Van der Waal broadening')
# plt.xlabel('Height [km]')
# plt.ylabel('$\gamma_{vdW}$ [s$^{-1}$]')
# plt.show()


# Compute damping parameter
wavelen = np.linspace(wav[0]-2e-8, wav[0]+2e-8, 1001)  # range around center

# prepare arrays
a = np.zeros((len(height), len(wavelen)))
v = np.zeros((len(height), len(wavelen)))
voigtprofile = np.zeros((len(height), len(wavelen)))

# calculate the voigt profile for the different height
gamma = gammavdw_nad(temp, pgas, 2, wav, E_ionization)  # van der waal damping
for j in range(len(height)):
    for i in range(len(wavelen)):
        a[j][i] = wavelen[i]**2/(4*np.pi*c)*gamma[j]/doppler[0][j]
        v[j][i] = (wavelen[i]-wav[0])/doppler[0][j]
        voigtprofile[j][i] = voigt(a[j][i], v[j][i])*doppler[0][j]*np.sqrt(np.pi)


deltawav = (wavelen - wav[0])*1e8  # delta wav in angstrom

# plt.figure(42)
# plt.plot(deltawav, voigtprofile[np.where(height == 0)][0], '-', label='0km')
# plt.plot(deltawav, voigtprofile[np.where(height == 200)][0], '--', label='200km')
# plt.plot(deltawav, voigtprofile[np.where(height == 400)][0], '-.', label='400km')
# plt.yscale('log')
# plt.xlabel('$\Delta \lambda$[$\AA$]')
# plt.title('Voigt function for Na I in FALC')
# plt.xlim([np.min(deltawav), np.max(deltawav)])
# plt.legend()
# plt.show()

correction = 1 - b_u/b_l*np.exp(-h*c/wav[0]/kerg/temp)

plt.figure(5)
plt.plot(height, correction)
plt.xlim([np.min(height), 2000])
plt.title(r'Correction term')
plt.xlabel(r'Height [km]')
plt.ylabel(r'')
plt.show()

# NaID1 extiction

NaD1_extinction = np.zeros((len(wavelen), len(height)))
ionstage = 1
level = 2
exthminD1 = np.zeros((len(wavelen), len(height)))
extcont = np.zeros((len(wavelen), len(height)))
exttotal = np.zeros((len(wavelen), len(height)))
thomson = 6.648*1e-25  # cm**2
extincelec = thomson*nel
for j in range(len(wavelen)):
    for i in range(len(height)):
        # Calculate NaID1 extinction
        part1 = np.sqrt(np.pi)*np.exp(2)/m_e/c*wavelen[j]*wavelen[i]/c*b_l
        part2 = sahabolt_na(temp[i], nel[i], ionstage, level, E_ionization)
        part3 = nhyd[i]*A_Na*f_lu[0]
        part4 = np.sqrt(np.pi)*voigtprofile[i, j]
        part5 = 1. - b_u/b_l*np.exp(-h*c/wav[0]/kerg/temp[i])
        NaD1_extinction[j][i] = part1*part2*part3*part4*part5
        # Calculate continuum extinction
        exthminD1[j][i] = exthmin(wavelen[j]*1e8, temp[i], nel[i])*(nhyd[i] - nprot[i])
        extcont[j][i] = exthminD1[j][i]+extincelec[i]
        # Sum to get the total
        exttotal[j][i] = NaD1_extinction[j][i] + extcont[j][i]

# plt.figure(43)
# plt.plot(height, NaD1_extinction[np.where(wavelen == wav[0])][0], '-', label='NaID1')
# plt.plot(height, extcont[np.where(wavelen == wav[0])][0], '--', label='Continuum')
# # plt.plot(height, exttotal[np.where(wavelen == wav[0])][0], '-.', label='Sum')
# plt.xlim([np.min(height), 2000])
# plt.title(r'Extinction at line center')
# plt.xlabel(r'Height [km]')
# plt.ylabel(r'Extinction $\alpha_\lambda^l$ [cm$^{-1}$]')
# plt.yscale('log')
# plt.legend(loc='best')
# plt.ylim([1e-14, 1e-2])
# plt.show()


# 3.7 Computed Na D1 line profile

# Calculate extiction for each wavelength at each height.
# Sum the one from the line extinction and the one from the continuum.
# Plot this against wavelength for selected heights

NaD1_extinction_swapped = np.zeros((len(height), len(deltawav)))
extcont_swapped = np.zeros((len(height), len(deltawav)))
exttotal_swapped = np.zeros((len(height), len(deltawav)))
for i in range(len(height)):
    for j in range(len(deltawav)):
        NaD1_extinction_swapped[i][j] = NaD1_extinction[j][i]
        extcont_swapped[i][j] = extcont[j][i]
        exttotal_swapped[i][j] = exttotal[j][i]

# plt.figure(44)
# plt.plot(deltawav, NaD1_extinction_swapped[np.where(height == 0)][0], '-.', label='NaI h=0')
# plt.plot(deltawav, extcont_swapped[np.where(height == 0)][0], '--', label='Cont h=0')
# plt.plot(deltawav, exttotal_swapped[np.where(height == 0)][0], label=r'Total h=0')
# plt.plot(deltawav, NaD1_extinction_swapped[np.where(height == 560)][0], '-.', label='NaI h=560')
# plt.plot(deltawav, extcont_swapped[np.where(height == 560)][0], '--', label='Cont h=560')
# plt.plot(deltawav, exttotal_swapped[np.where(height == 560)][0], label=r'Total h=560')
# plt.yscale('log')
# plt.xlim([np.min(deltawav), np.max(deltawav)])
# plt.title(r'Extinction')
# plt.ylabel(r'Extinction $\alpha_\lambda^l$ [cm$^{-1}$]')
# plt.xlabel(r'$\Delta\lambda$[$\AA$]')
# plt.legend(loc='best')
# plt.show()


# Compute the extinction for NaID1line
tau = np.zeros(np.size(tau5))
integrand = np.zeros(np.size(tau5)) 
contfunc = np.zeros(np.size(tau5))

intt = np.zeros(np.size(wavelen))
hint = np.zeros(np.size(wavelen))

for j in range(np.size(wavelen)):
    for i in range(1, len(tau)):  # the index zero is not accounted for
        tau[i] = tau[i-1] + 0.5*(exttotal[j][i]+exttotal[j][i-1])*(height[i-1]-height[i])*1e5
        integrand[i] = planck(temp[i], wavelen[j]*1e4)*np.exp(-tau[i])
        intt[j] += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
        hint[j] += height[i]*0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])

wavelen *= 1e8  # convert wavelen from cm to Angstrom

plt.figure(45)
plt.plot(wavelen, intt/np.max(intt), '-')
plt.plot(wavelength_air, S_spec, '--')
plt.xlim([wav[0]*1e8-2, wav[0]*1e8+2])
plt.title(r'Disk center intensity around line center of the NaID1 line')
plt.xlabel(r'Wavelength[$\AA$]')
plt.ylabel(r'Normalized intensity')
ax = plt.gca()  # This disables the offset
ax.ticklabel_format(useOffset=False)  # disable offset
plt.show()
