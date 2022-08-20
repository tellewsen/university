"""
Term Project 2
The second term project now involves using Chap.5 to include convection
in your model of a star. In the first term project you modeled the sun from
its core to the edge if its radiative zone. Now you can include the rest of the
sun by including energy transport through convection. You should use the
model you have already done, but now include a check at each radius (or
mass) shell of convective stability. In case the shell is not convectively stable,
you need to use the equation for convective transport instead of radiative
transport. That involves solving Exerc. 5.10, to get the expression for the
convective flux, and through that get the temperature gradient. You can
make the same assumptions as in the first term project, for instance only
include PPI and PPII etc.

The project can be found at 
http://www.uio.no/studier/emner/matnat/astro/AST3310/v15/notes/exerc2.pdf

Created by Andreas Ellewsen
"""
import numpy as np  # Numerical python
from matplotlib.pyplot import plot, figure, title, xlabel, ylabel, show, legend, savefig  # Plotting library
from functions import energy, density, pressure, opacity, xi_func

# Constants
sigma = 5.67E-8             # [W*m**-2*K**-4] Stefan-Boltzmann's constant
G = 6.67384E-11	            # [m**3*kg**-1*s**-2] Gravitational constant
MeVtoJ = 1.60217657E-13     # [J] Conversion factor from MeV to Joules
L_sun = 3.846E26            # [W] Luminosity of the sun
R_sun = 6.958E8	            # [m] Radius of the sun
M_sun = 1.989E30            # [kg] Mass of the sun
alpha = 1

# Initial conditions
X = .7              # [] Hydrogen fraction
Y3 = 1E-10          # [] Helium 3 fraction
Y = .29             # [] Sum of Helium 3 and Helium 4 fractions
Y4 = Y-Y3           # [] Helium 4 fraction
Z_73Li = 1E-13      # [] Lithium  fraction
Z_74Be = 1E-13      # [] Beryllium fraction
Z = .01		        # [] Other
L_0 = L_sun         # [W] Luminosity of the star
R_0 = R_sun   	    # [m] Radius of the star
M_0 = M_sun         # [kg] Mass of star
rho_0 = 4E-4        # [kg*m**-3] Mass density of star
P_0 = 1E8           # [Pa] Test value for pressure
T_0 = 5770	  	    # [K] Temperature
nabla_ad = 2./5     # Adiabatic gradient for fully ionized ideal gas
mu = 1. / (2. * X + 7. / 4 * Y + 5. / 7 * Z_74Be + 4. / 7 * Z_73Li + Z)


# Solving the equations dictating the physics of the sun.
'''This is the part of the program that actually solves the equations\
governing the sun.'''

# Make empty lists
R = []
L = []
P = []
T = []
M = []
rho = []
epsilon = []
F_C = []
F_R = []
Convection_frac = []
Radiation_frac = []
PPI_fraction = []
PPII_fraction = []
PPI_only = []
PPII_only = []
PPI_frac = []
PPII_frac = []

# Fill the lists with the initial conditions
M.append(M_0)
R.append(R_0)
L.append(L_0)
T.append(T_0)
rho.append(rho_0)
P.append(pressure(rho[0], T[0], mu))
epsilon.append(energy(T[0], rho[0], X, Y3, Y4, Z_73Li, Z_74Be)[0])
# Simulate evolution by reducing radius.
i = 0  # Starts a counter
p = 0.01  # Max change allowed per step.

while M[i] > 0 and R[i] > 0 and L[i] > 0 and T[i] > 0 and P[i] > 0:
    # Change in radius
    dMdr = 4*np.pi*R[i]**2*rho[i]
    drM = p*M[i]/dMdr   # Largest dr allowed

    # Change in Pressure
    dPdr = -(G*M[i]*rho[i])/R[i]**2
    drP = p*P[i]/dPdr    # Largest dr allowed

    # Change in Luminosity
    dLdr = epsilon[i]*dMdr
    drL = p*L[i]/dLdr  # Largest dr allowed

    # Calculate the gradients
    lines = np.genfromtxt('opacity.txt')
    kappa = opacity(T[i], rho[i], lines)  # Opacity

    g = G*M[i]/R[i]**2  # Gravitational acceleration at radius R
    c_P = 5./2*P[i]/(rho[i]*T[i])  # Heat capacity constant pressure
    H_P = P[i]*R[i]**2 / (G*M[i]*rho[i])  # Pressure scale height
    l_m = alpha*H_P  # Mixing length
    U = 64*sigma*T[i]**3/(3*kappa*rho[i]**2*c_P)*(H_P/g)**(1./2)
    nabla_rad = 3*kappa*L[i]*P[i]/(64*np.pi*sigma*G*T[i]**4*M[i])
    
    # Calculate xi
    xi = xi_func(l_m, U, nabla_ad, nabla_rad)
    # Calculate variables depending on xi
    # v = (g*l_m**2/(4*H_P))**(1./2)*xi
    nabla = nabla_rad - l_m**2/U*xi**3
    nabla_star = nabla - xi**2

    """
    #Print test values
    print 'nabla_rad  = ',nabla_rad
    print 'nabla_ad   = ',nabla_ad
    print 'nabla      = ',nabla
    print 'nabla_star = ',nabla_star
    """
    
    # Change in temperature
    'Instability criterion deciding which gradient to use.'
    if nabla_rad > nabla_ad:
        dTdr = -T[i]/H_P*nabla
        F_R.append(-16./3*sigma*T[i]**3/(kappa*rho[i])*dTdr)
        F_C.append(L[i]/(4*np.pi*R[i]**2)-F_R[i])
    else:
        dTdr = -T[i]/H_P*nabla_rad
        F_C.append(0)
        F_R.append(-16./3*sigma*T[i]**3/(kappa*rho[i])*dTdr)
    drT = p*T[i]/dTdr  # Largest dr allowed
    
    # Decide how small dr needs to be
    drR = 10*p*R[i]
    dr = -min(abs(drR), abs(drP), abs(drL), abs(drT), abs(drM))
    if abs(dr) < 1E2:
        dr = -1E2   # Limit how small dr is allowed to become

    # Update values
    M.append(M[i] + dMdr*dr)		# Update Mass
    R.append(R[i] + dr)     		# Update Radius
    L.append(L[i] + dLdr*dr) 		# Update Luminosity
    P.append(P[i] + dPdr*dr)		# Update Pressure
    T.append(T[i] + dTdr*dr)		# Update Temperature
    rho.append(density(P[i], T[i], mu))	 # Update Density
    epsilon.append(energy(T[i], rho[i], mu, Y3, Y4, Z_73Li, Z_74Be)[0])  # Update Energy output
    PPI_fraction.append(energy(T[i], rho[i], mu, Y3, Y4, Z_73Li, Z_74Be)[1])  # PPI of total
    PPII_fraction.append(energy(T[i], rho[i], mu, Y3, Y4, Z_73Li, Z_74Be)[2])  # PPII of total
    PPI_only.append(energy(T[i], rho[i], mu, Y3, Y4, Z_73Li, Z_74Be)[3])  # PPI without start
    PPII_only.append(energy(T[i], rho[i], mu, Y3, Y4, Z_73Li, Z_74Be)[4])  # PPII without start
    Convection_frac.append(F_C[i]/(F_C[i]+F_R[i]))  # Convection fraction
    Radiation_frac.append(F_R[i]/(F_C[i]+F_R[i]))  # Radiation fraction
    
    """print 'M = %.3e, R = %.3e, L = %.3e, P = %.3e, T= %.3e, rho = %.3e, 
    kappa = %.3e, epsilon = %.3e, dT = %.3e, dr = %.3e, dL = %.3e, dm = %.3e, 
    dP = %.3e, F_C = %.3e, F_R = %.3e'%(M[i], R[i], L[i], P[i], T[i], rho[i], 
    kappa, epsilon[i], dTdr*dr, dr, dLdr*dr, dMdr*dr, dPdr*dr,F_C,F_R)"""

    # Print  information about which variable went to zero first.
    if M[i+1] < 0:
        print 'Mass went negative'
    if R[i+1] < 0:
        print 'Radius went negative'
    if L[i+1] < 0:
        print 'Luminosity went negative'
    if T[i+1] < 0:
        print 'Temperature went negative'
    if P[i+1] < 0:
        print 'Pressure went negative'

    i += 1   # increase counter

# Fractions of the total production for each of the chains.
PPI_fraction = np.array(PPI_fraction)
PPII_fraction = np.array(PPII_fraction)
PPI_frac = PPI_fraction/(PPI_fraction + PPII_fraction)
PPII_frac = PPII_fraction/(PPI_fraction + PPII_fraction)

# Fraction of the total production for each of the chains without the prerequisites.
PPI_only = np.array(PPI_only)
PPII_only = np.array(PPII_only)
PPI_only_frac = PPI_only/(PPI_only + PPII_only)
PPII_only_frac = PPII_only/(PPI_only + PPII_only)

# Removes the lat set of results since they are unphysical
del M[-1]
del R[-1]
del L[-1]
del P[-1]
del T[-1]
del rho[-1]
del epsilon[-1]

# Print start values used in this run
print 'Start values:'
print 'M_0 :		   %.3e' % (M[0])
print 'R_0 :		   %.3e' % (R[0])
print 'L_0 :		   %.3e' % (L[0])
print 'P_0 :		   %.3e' % (P[0])
print 'T_0 :		   %.3e' % (T[0])
print 'rho_0 :		   %.3e' % (rho[0])
print 'epsilon_0:	   %.3e' % (epsilon[0])
print 'F_C: ', F_C[0]
print 'F_R: ', F_R[0]

# Scale values after initial conditions
M = np.array(M)/M[0]
R = np.array(R)/R[0]
L = np.array(L)/L[0]

# print Final values
print 'Final values:'
print 'M/M_0 :		   %.3f' % (M[-1])
print 'R/R_0 :		   %.3f' % (R[-1])
print 'L/L_0 :		   %.3f' % (L[-1])
print 'P : 		   %.3e' % (P[-1])
print 'T :  		   %.3e' % (T[-1])
print 'rho   :       %.3e' % (rho[-1])
print 'epsilon:	   %.3e' % (epsilon[-1])
print 'F_C: ', F_C[-1]
print 'F_R: ', F_R[-1]

# Plots of interesting things
figure(0)
plot(R, M)
title('Mass vs Radius')
xlabel('Radius/R_0')
ylabel('Mass[kg]/M_0')
savefig('figures/mass.png')

figure(1)
plot(R, P)
title('Pressure vs Radius')
xlabel('Radius/R_0')
ylabel('Pressure[Pa]')
savefig('figures/pressure.png')

figure(2)
plot(R, L)
title('Luminosity vs Radius')
xlabel('Radius/R_0')
ylabel('Luminosity/L_0')
savefig('figures/luminosity.png')

figure(3)
plot(R, T)
title('Temperature vs Radius')
xlabel('Radius/R_0')
ylabel('Temperature[K]')
savefig('figures/temperature.png')

figure(4)
plot(R, rho)
title('Density vs Radius')
xlabel('Radius/R_0')
ylabel('Density[kg/m^3]')
savefig('figures/density.png')

figure(5)
plot(R, Convection_frac)
title('Energy transport vs Radius')
xlabel('Radius/R_0')
ylabel('Energy/Total Energy')
plot(R, Radiation_frac)
legend(['Convection', 'Radiation'])
savefig('figures/energytransport.png')

figure(6)
plot(R, epsilon)
title('Energy production vs Radius')
xlabel('Radius/R_0')
ylabel('Energy[J/m]')
plot(R, PPI_fraction)
plot(R, PPII_fraction)
legend(['Total', 'PPI chain', 'PPII Chain'])
savefig('figures/energyproduction.png')

figure(7)
plot(R, F_C)
title('Energy transport vs Radius')
xlabel('Radius/R_0')
ylabel('Energy[W/m^2]')
plot(R, F_R)
legend(['Convection', 'Radiation'])
savefig('figures/energytransport_frac.png')

figure(8)
plot(R, PPI_frac)
title('Energy production vs Radius')
xlabel('Radius/R_0')
ylabel('Energy/Total energy')
plot(R, PPII_frac)
legend(['PPI', 'PPII'])
savefig('figures/energyproduction_frac.png')

figure(9)
plot(R, PPI_only_frac)
title('Energy production vs Radius')
xlabel('Radius/R_0')
ylabel('Energy/Total energy')
plot(R, PPII_only_frac)
legend(['PPI', 'PPII'])
savefig('figures/energyproduction_frac_only.png')

show()
