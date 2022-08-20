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
from numpy import exp,zeros\
,genfromtxt,log10,pi,array,roots,imag,real    #Numerical python

from matplotlib.pyplot import plot,figure\
,title,xlabel,ylabel,show,hold,legend,savefig	#Plotting library

import sys as sys	#System functions

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
alpha  = 1
#Energy output of reactions
Q_pp = (.15+1.02)  *MeVtoJ	#[J]
Q_dp = (5.49)	   *MeVtoJ	#[J]
Q_33 = (12.86)	   *MeVtoJ	#[J]
Q_34 = (1.59)	   *MeVtoJ	#[J]
Q_7e = (.05)	   *MeVtoJ	#[J]
Q_71prime = (17.35)*MeVtoJ	#[J]

#Initial conditions
X       = .7      #[] Hydrogen fraction
Y3      = 1E-10   #[] Helium 3 fraction
Y       = .29     #[] Sum of Helium 3 and Helium 4 fractions
Y4      = Y-Y3    #[] Helium 4 fraction
Z_73Li  = 1E-13   #[] Lithium  fraction
Z_74Be  = 1E-13   #[] Berylium fraction
Z       = .01		#[] Other
L_0     = L_sun   #[W] Luminosity of the star
R_0     = R_sun   	#[m] Radius of the star
M_0     = M_sun	#[kg] Mass of star 
rho_0   = 4E-4		#[kg*m**-3] Mass density of star 
P_0     = 1E8     #[Pa] Test value for pressure
T_0     = 5770	  	#[K] Temperature
nabla_ad = 2./5   #Adiabatic gradient for fully ionized ideal gas

#Computes number of atoms of the different elements.
def atomnumbers(rho):
	n_e   = rho*(1+X) /(2.*m_u)    	#Number of electrons
	n_d   = 0			#Number of deuterium atoms
	n_p   = rho*X	  /m_u     	#Number of Hydrogen atoms
	n_3   = rho*Y3	  /(3.*m_u)     #Number of Helium 3 atoms
	n_4   = rho*Y4	  /(4.*m_u)     #Number of Helium 4 atoms
	n_Li  = rho*Z_73Li/(7.*m_u) 	#Numver of Beryllium 7 atoms
	n_Be  = rho*Z_74Be/(7.*m_u) 	#Number of Lithium 7 atoms
	atomnumbers = [n_e,n_d,n_p,n_3,n_4,n_Li,n_Be]
	return atomnumbers

#Reactions
"""In this part the reaction rates of the processes in the PPI and PPII \
chain are calculated. I neglect the other processes since I choose to \
look at a star with low temperature at it's core, and in that case those \
processes are very slow compared to PPI and PPII."""

def reactions(T):
	T9 		= T*1E-9		#Convert T to form used in reactions.
	T9_star1 	= T9/(1+4.95E-2*T9)	#Other form used.
	T9_star2 	= T9/(1+.759*T9)	#Yet another form used.

	H_D   	= 	4.01E-15*T9**(-2./3)*exp(-3.380*T9**(-1./3))*(1 +\
                     .123*T9**(1./3) + 1.09*T9**(2./3) + .938*T9)

	He3_p	= 	6.04E10*T9**(-2./3)*exp(-12.276*T9**(-1./3))*(1 +\
                     .034*T9**(1./3) - .522*T9**(2./3)-.124*T9 + \
                     .353*T9**(4./3)+ .213*T9**(-5./3))

	He3_Be  = 	5.61E6*T9_star1**(5./6)*T9**(-3./2)*exp(-12.826*\
                     T9_star1**(-1./3))

	Be_Li   = 	1.34E-10*T9**(-1./2)*(1 - .537*T9**(1./3)  + \
                     3.86*T9**(2./3) +  .0027*T9**-1*exp(2.515E-3*T9**-1))

	Li_p	= 	1.096E9*T9**(-2./3)*exp(-8.472*T9**(-1./3)) - \
                     4.830E8*T9_star2**(5./6)*T9**(-3./2)*exp(-8.472*\
                     T9_star2**(-1./3)) + 1.06E10*T9**(-3./2)\
                     *exp(-30.442*T9**-1)
	"""
	"These two reactions are not used in this assignment, but may be useful \
     at a later point, so I've kept them here."
     
	Berylliym_Boron = 	3.11E5*T9**(-2./3)*exp(-10.262*T9**(-1./3))+2.53E3\
                         *T9**(-3./2)*exp(-7.306*T9**-1)

	Nitrogen_Oxygen = 	4.9E7*T9**(-2./3)*exp(-15.228*T9**(-1./3)-\
                         .092*T9**2)*(1+.027*T9**(1./3)-.778*T9**(2./3)\
                         -.149*T9+.261*T9**(4./3)+.127*T9**(5./3)) + \
                         2.37E3*T9**(2./3)*exp(-3.011*T9**-1) + \
                         2.19E4*exp(-12.53*T9**-1)
	"""
	rrs = [H_D,He3_p,He3_Be,Be_Li,Li_p]#,Berylliym_Boron,Nitrogen_Oxygen
	return rrs

#Computing reactionrates
"""
This segment computes the reactionrates for the different processes in the sun.
"""
def Reactionrates(T,rho):
	n = atomnumbers(rho)
	rrs = reactions(T)
 
	#Reaction rates
	"""All of the following rates are calculated by the equation
	r_ik = n_i*n_k/(rho*(1+delta_ik)*lambda_ik)
	where i,k define the elements and delta_ik =1 if i=k else 0.
	"""
	lambda_pp = rrs[0]/N_A*1E-6		#[m^3/s]
	lambda_33 = rrs[1]/N_A*1E-6		#[m^3/s]
	lambda_34 = rrs[2]/N_A*1E-6		#[m^3/s]
	lambda_7e = rrs[3]/N_A*1E-6		#[m^3/s]
	lambda_71prime  = rrs[4]/N_A*1E-6	#[m^3/s]

	r_pp 	  = n[2]**2/(rho*2)*lambda_pp	#[kg-1*s-1]
	r_pd  	  = r_pp			#[kg-1*s-1] Assume this reaction happens \
						#instantly such that it happens at the same \
						#rate as the elements it needs become available.
						#Thus it must be the same as the reaction that \
						#creates those elements (in this case Deuterium)
	r_33 	  = n[3]**2  /(rho*2)*lambda_33	#[kg-1*s-1]
	r_34 	  = n[3]*n[4]/rho*lambda_34	#[kg-1*s-1]
	r_7e 	  = n[5]*n[6]/rho*lambda_7e	#[kg-1*s-1]
	r_71prime = n[5]*n[2]/rho*lambda_71prime#[kg-1*s-1]
		
	#This part makes sure that reaction rates which rely on other rates \
	#don't use elements that are not yet present.
	if 2.*r_33 + r_34 > r_pd:
		r_33 = 2./3*r_pd
		r_34 = 1./3*r_pd
	if r_7e > r_34:
		r_7e = r_34
	if r_71prime > r_7e:
		r_71prime = r_7e
	reactionrate = r_pp,r_pd,r_33,r_34,r_7e,r_71prime
	return reactionrate

#Solving the equation for the change in Luminosity

def Energy(T,rho):
	reactionrate = Reactionrates(T,rho)
	e_1  = reactionrate[0]*(Q_pp+Q_dp)
	e_2  = reactionrate[2]*Q_33
	e_3  = reactionrate[3]*Q_34
	e_4  = reactionrate[4]*Q_7e
	e_5  = reactionrate[5]*Q_71prime
	e    = e_1 + e_2 + e_3 + e_4 + e_5
	PPI_frac  = e_2 + .69*e_1
	PPII_frac = e_3 + e_4 + e_5 + .31*e_1 
	PPI_only = e_2
	PPII_only = e_3 +e_4+e_5
	output = e,PPI_frac,PPII_frac,PPI_only,PPII_only
	return output

#Defining functions for pressure,temperature and density
"""This part solves the equation for Pressure, assuming 
the equation of state to be that of an ideal gas"""

#Constants needed
mu 	= 1/(2.*X + 7./4*Y + 5./7*Z_74Be + 4./7*Z_73Li + Z)
a 	= 4.*sigma/c

def Pressure(rho,T):
	P_rad 	= a/.3*T**4           #Radiation pressure
	P_G 	= rho*k*T/(mu*m_u)    #Gas pressure, assuming ideal gas
	P	= P_G + P_rad            #Total Pressure
	return P

def Density(P,T):
	rho  = (P*mu*m_u)/(k*T)		#Ideal gas
	return rho

#Function that reads the opacity table
lines = genfromtxt('opacity.txt')

def Opacity(T,rho):
	'''This function reads the opacity table'''
	#Convert R and T to the form used in opacity.txt
	R      = (rho*1E-3)/(T/1E6) 
	Rvalue = log10(R)
	Tvalue = log10(T)
	Rvalue = (round(2*Rvalue)/2.)
	#Print errormessage if using values outside table
	if Rvalue >= 1.5 or Rvalue <= -8.5:
		print 'Tried using R outside opacity table'
		sys.exit()
	if Tvalue >= 8.8 or Tvalue <= 3.7 :
		print 'Tried using T outside opacity table'
		sys.exit()

	#Pick the right R value from the table
	for i in range(len(lines[0])):
		if Rvalue == lines[0][i]:
			Ri = i

   	#Make list of T values found in table
	Tlist = zeros(len(lines))
	for i in range(1,len(lines)):
		Tlist[i] = lines[i][0]

	#Pick the right T value from the table
	Ti = min(range(len(Tlist)), key=lambda i: abs(Tlist[i]-Tvalue))

	#Pick the right kappa value with R and T
	log10kappa = lines[Ti][Ri]

	#Convert the value given to the form used in the equations
	kappa = .1*10.**(log10kappa)
	return kappa

#Solving the equations dictating the physics of the sun.
'''This is the part of the program that actually solves the equations\
governing the sun.'''

#Make empty lists
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
PPII_only =[]
PPI_frac = []
PPII_frac = []

#Fill the lists with the initial conditions
M.append(M_0)
R.append(R_0)
L.append(L_0)
T.append(T_0)
rho.append(rho_0)
P.append(Pressure(rho[0],T[0]))
epsilon.append(Energy(T[0],rho[0])[0])


#Function to calculate xi
def xifunc(l_m,U,nabla_ad,nabla_rad):
    c1 = l_m**2 / U
    c2 = 1
    c3 = 4*U/l_m**2
    c4 = nabla_ad -nabla_rad
    xi = roots([c1,c2,c3,c4])
    for element in xi:
        if imag(element) == 0:
            xi = real(element)
            break
    return xi
    
#Simulate evolution by reducing radius.
i = 0    #Starts a counter
p = 0.01 #Max change allowed per step.

while M[i] > 0 and R[i] > 0 and L[i] > 0 and T[i] > 0 and P[i] > 0:
    #Change in radius
    dMdr	= 4*pi*R[i]**2*rho[i]
    drM 	= p*M[i]/dMdr   #Largest dr allowed

    #change in Pressure
    dPdr 	= -(G*M[i]*rho[i])/R[i]**2
    drP 	= p*P[i]/dPdr    #Largest dr allowed

    #Change in Luminosity
    dLdr 	= epsilon[i]*dMdr
    drL 	= p*L[i]/dLdr  #Largest dr allowed

    #Calculate the gradients
    g       = G*M[i]/R[i]**2    #Gravitational acceleration at radius R
    kappa   = Opacity(T[i],rho[i])		#Opacity
    c_P     = 5./2*P[i]/(rho[i]*T[i])	#Heat capacity constant pressure
    H_P     = P[i]*R[i]**2 / (G*M[i]*rho[i])	#Pressure scale height
    l_m     = alpha*H_P	                #Mixing length
    U       = 64*sigma*T[i]**3/(3*kappa*rho[i]**2*c_P)*(H_P/g)**(1./2)
    nabla_rad  = 3*kappa*L[i]*P[i]/(64*pi*sigma*G*T[i]**4*M[i])
    
    #Calculate xi
    xi = xifunc(l_m,U,nabla_ad,nabla_rad)
    #Calculate variables depening on xi
    #v = (g*l_m**2/(4*H_P))**(1./2)*xi
    nabla      = nabla_rad - l_m**2 / U *xi**3
    nabla_star = nabla - xi**2


    """#Print test values
    print 'nabla_rad  = ',nabla_rad
    print 'nabla_ad   = ',nabla_ad
    print 'nabla      = ',nabla
    print 'nabla_star = ',nabla_star"""
    
    #Change in temperature    
    'Instability criterion deciding which gradient to use.'
    if nabla_rad > nabla_ad:
        dTdr = -T[i]/H_P*nabla
        F_R.append(-16./3*sigma*T[i]**3/(kappa*rho[i])*dTdr)
        F_C.append(L[i]/(4*pi*R[i]**2)-F_R[i])
    else:
        dTdr = -T[i]/H_P*nabla_rad
        F_C.append(0)
        F_R.append(-16./3*sigma*T[i]**3/(kappa*rho[i])*dTdr)
    drT 	= p*T[i]/dTdr  #Largest dr allowed
    
    #Decide how small dr needs to be
    drR =   10*p*R[i]
    dr  =   -min(abs(drR),abs(drP),abs(drL),abs(drT),abs(drM))
    if abs(dr) < 1E2:
        dr = -1E2   #Limit how small dr is allowed to become

   #Update values
    M.append(M[i] + dMdr*dr)		#Update Mass
    R.append(R[i] + dr)     		#Update Radius
    L.append(L[i] + dLdr*dr) 		#Update Luminosity
    P.append(P[i] + dPdr*dr)		#Update Pressure
    T.append(T[i] + dTdr*dr)		#Update Temperature
    rho.append(Density(P[i],T[i]))	#Update Density
    epsilon.append(Energy(T[i],rho[i])[0]) #Update Energy output
    PPI_fraction.append(Energy(T[i],rho[i])[1]) #PPI of total
    PPII_fraction.append(Energy(T[i],rho[i])[2]) #PPII of total
    PPI_only.append(Energy(T[i],rho[i])[3]) #PPI without start
    PPII_only.append(Energy(T[i],rho[i])[4]) #PPII without start
    Convection_frac.append(F_C[i]/(F_C[i]+F_R[i])) #Convection fraction
    Radiation_frac.append(F_R[i]/(F_C[i]+F_R[i])) #Radiation fraction
    
    """print 'M = %.3e, R = %.3e, L = %.3e, P = %.3e, T= %.3e, rho = %.3e, \
kappa = %.3e, epsilon = %.3e, dT = %.3e, dr = %.3e, dL = %.3e, dm = %.3e, \
dP = %.3e, F_C = %.3e, F_R = %.3e'%(M[i], R[i], L[i], P[i], T[i], rho[i], \
kappa, epsilon[i], dTdr*dr, dr, dLdr*dr, dMdr*dr, dPdr*dr,F_C,F_R)"""

    #Print  information about which variable went to zero first.   
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

    i +=1   #increase counter

#Fractions of the total production for each of the chains.
PPI_fraction = array(PPI_fraction)
PPII_fraction = array(PPII_fraction)
PPI_frac = PPI_fraction/(PPI_fraction + PPII_fraction)
PPII_frac = PPII_fraction/(PPI_fraction + PPII_fraction)

#Fraction of the total prodcution for each of the chains without
#the prerequisites.
PPI_only = array(PPI_only)
PPII_only = array(PPII_only)
PPI_only_frac = PPI_only/(PPI_only + PPII_only)
PPII_only_frac = PPII_only/(PPI_only + PPII_only)

#Removes the lat set of results since they are unphysical
del M[-1]
del R[-1]
del L[-1]
del P[-1]
del T[-1]
del rho[-1]
del epsilon[-1]

#Print start values used in this run
print 'Start values:'
print 'M_0 :		   %.3e'%(M[0])
print 'R_0 :		   %.3e'%(R[0])
print 'L_0 :		   %.3e'%(L[0])
print 'P_0 :		   %.3e'%(P[0])
print 'T_0 :		   %.3e'%(T[0])
print 'rho_0 :		   %.3e'%(rho[0])
print 'epsilon_0:	   %.3e'%(epsilon[0])
print 'F_C: ', F_C[0]
print 'F_R: ', F_R[0]

#Scale values after initial conditions
M = array(M)/M[0]
R = array(R)/R[0]
L = array(L)/L[0]

#print Final values
print 'Final values:'
print 'M/M_0 :		   %.3f'%(M[-1])
print 'R/R_0 :		   %.3f'%(R[-1])
print 'L/L_0 :		   %.3f'%(L[-1])
print 'P : 		   %.3e'%(P[-1])
print 'T :  		   %.3e'%(T[-1])
print 'rho   :       %.3e'%(rho[-1])
print 'epsilon:	   %.3e'%(epsilon[-1])
print 'F_C: ', F_C[-1]
print 'F_R: ', F_R[-1]

#Plots of interesting things
figure()
plot(R,M)
title('Mass vs Radius')
xlabel('Radius/R_0')
ylabel('Mass[kg]/M_0')
savefig('mass.png')

figure()
plot(R,P)
title('Pressure vs Radius')
xlabel('Radius/R_0')
ylabel('Pressure[Pa]')
savefig('pressure.png')

figure()
plot(R,L)
title('Luminosity vs Radius')
xlabel('Radius/R_0')
ylabel('Luminosity/L_0')
savefig('luminosity.png')

figure()
plot(R,T)
title('Temperature vs Radius')
xlabel('Radius/R_0')
ylabel('Temperature[K]')
savefig('temperature.png')

figure()
plot(R,rho)
title('Density vs Radius')
xlabel('Radius/R_0')
ylabel('Density[kg/m^3]')
savefig('density.png')

figure()
plot(R,Convection_frac)
title('Energy transport vs Radius')
xlabel('Radius/R_0')
ylabel('Energy/Total Energy')
hold('on')
plot(R,Radiation_frac)
legend(['Convection','Radiation'])
hold('off')
savefig('energytransport.png')

figure()
plot(R,epsilon)
title('Energy production vs Radius')
xlabel('Radius/R_0')
ylabel('Energy[J/m]')
hold('on')
plot(R,PPI_fraction)
plot(R,PPII_fraction)
legend(['Total','PPI chain','PPII Chain'])
hold('off')
savefig('energyproduction.png')

figure()
plot(R,F_C)
title('Energy transport vs Radius')
xlabel('Radius/R_0')
ylabel('Energy[W/m^2]')
hold('on')
plot(R,F_R)
legend(['Convection','Radiation'])
hold('off')
savefig('energytransport_frac.png')

figure()
plot(R,PPI_frac)
title('Energy production vs Radius')
xlabel('Radius/R_0')
ylabel('Energy/Total energy')
hold('on')
plot(R,PPII_frac)
legend(['PPI','PPII'])
hold('off')
savefig('energyproduction_frac.png')


figure()
plot(R,PPI_only_frac)
title('Energy production vs Radius')
xlabel('Radius/R_0')
ylabel('Energy/Total energy')
hold('on')
plot(R,PPII_only_frac)
legend(['PPI','PPII'])
hold('off')
savefig('energyproduction_frac_only.png')

show()