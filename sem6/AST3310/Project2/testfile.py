from numpy import pi,roots,real
#Constants
G	= 6.67384E-11		#[m**3*kg**-1*s**-2] Gravitational constant
m_u 	= 1.66053892173E-27 	#[kg] Atomic mass unit
N_A 	= 6.0221413E23      	#[] Avogadros number
k   	= 1.3806488E-23		#[m**2*kg*s**-2*K**-1] Boltzmans constant
c	= 2.99792458E8		#[m*s**-1] Speed of lights in vacuum
sigma   = 5.67E-8		#[W*m**-2*K**-4] Stefan-Boltzmans constant
MeVtoJ 	= 1.60217657E-13 	#[J] Conversion factor from MeV to Joules
L_sun  	= 3.846E26		#[W] Luminosity of the sun
R_sun  	= 6.958E8		#[m] Radius of the sun
M_sun  	= 1.989E30		#[kg] Mass of the sun

#init conditiosn
mu = .6
T = .9*10**6
rho  = 55.9
R    = .84*R_sun
M = .99*M_sun
kappa = 3.98
alpha = 1.
P =  rho*k*T/(mu*m_u)
L = L_sun
a 	= 4.*sigma/c

g = G*M/R**2
c_P = 5./2*P/(rho*T)	#Heat capacity constant pressure

H_P = P*R**2 / (G*M*rho)	#Pressure scale height
print 'H_P =', H_P

nabla_rad  = H_P/T*3.*kappa*L/(256.*pi**2.*sigma*R**4.*T**3.)*(4*pi*rho*R**2.)
print 'nabla_rad =' ,nabla_rad

nabla_ad = 2./5
print'nabla_ad =',nabla_ad

l_m  = alpha*H_P	#Mixing length

U = 64*sigma*T**3/(3*kappa*rho**2*c_P)*(H_P/g)**(1./2)
print 'U =',U

#Calculate xi
c1 = l_m**2 / U
c2 = 1
c3 = 4*U/l_m**2
c4 = nabla_ad -nabla_rad
xi = roots([c1,c2,c3,c4])

for i  in xi:
    if imag(i) == 0:
        xi = real(i)

print 'xi = ',xi


v = (g*l_m**2/(4*H_P))**(1./2)*xi
print 'v = ',v

F_C = rho*c_P*v*T*l_m/(2*H_P)*xi**2

nabla      = nabla_rad -3*kappa*P*R**2/(4*a*c*G*T**4*M)*F_C
print 'nabla = ',nabla

nabla_star = (64*sigma*T**3/(3*kappa*rho**2*c_P)*(H_P/g)**(1./2)*4/l_m**2 *xi) + nabla_ad
print 'nabla_star = ',nabla_star

F_rad = 4*a*c*G*T**4*M/(3*kappa*P*R**2)*nabla

Convection_frac = F_C/(F_C+F_rad)
print 'Convection_frac =',Convection_frac

Radiation_frac = F_rad/(F_C+F_rad)
print 'Radiation_frac =',Radiation_frac



