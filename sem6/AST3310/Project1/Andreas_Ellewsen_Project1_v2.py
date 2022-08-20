"""
Term Project 1
The first term project involves modeling the central parts of the sun by
writing a computer program to solve the governing equations. The project
involves information contained in Chaps.1-5.

The project can be found at 
http://www.uio.no/studier/emner/matnat/astro/AST3310/v15/notes/term1v1.pdf

Created by Andreas Ellewsen
"""
from numpy import exp,zeros,array_equal\
,genfromtxt,log10,pi    #Numerical python
from matplotlib.pyplot import plot,figure\
,title,xlabel,ylabel,show	#Plotting library
import sys as sys

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

#Energy output of reactions
Q_pp = (.15+1.02)  *MeVtoJ	#[J]
Q_dp = (5.49)	   *MeVtoJ	#[J]
Q_33 = (12.86)	   *MeVtoJ	#[J]
Q_34 = (1.59)	   *MeVtoJ	#[J]
Q_7e = (.05)	   *MeVtoJ	#[J]
Q_71prime = (17.35)*MeVtoJ	#[J]

#Initial conditions
X       = .7        	#[] Hydrogen fraction
Y3      = 1E-10     	#[] Helium 3 fraction
Y       = .29	    	#[] Sum of Helium 3 and Helium 4 fractions
Y4      = Y-Y3      	#[] Helium 4 fraction
Z_73Li  = 1E-13      	#[] Lithium  fraction
Z_74Be  = 1E-13     	#[] Berylium fraction
Z       = .01		#[] Other
L_0     =  L_sun    	#[W] Luminosity of the star
R_0     = .5*R_sun   	#[m] Radius of the star
M_0     = .7*M_sun	#[kg] Mass of star 
rho_0   = 1e3		#[kg*m**-3] Mass density of star 
T_0     = 1E5	  	#[K] Temperature
#P_0     = 1E11		#[Pa] Pressure


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

	Li_p	= 	     1.096E9*T9**(-2./3)*exp(-8.472*T9**(-1./3)) - \
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

#Time evolution
"""This segement takes care of the time evolution of the reactions. 
By doing this we hope to reach some kind of equilibrium over time, 
such that all the reaction rates stabilize at some level. 
This is achieved remarkably fast, as can be seen in the choice of Time 
and steps N."""
def Reactionrates(T,rho):
	counter = 0
	Time  = 1		#Total time to simulate [s]
	N  = 10**4		#Number of steps
	dt = Time/N		#Timesteps [s]
	n = atomnumbers(rho)
	rrs = reactions(T)
	reactionrate_old = zeros(6)
	reactionrate_new = zeros(6)
	for i in range(0,N-1):
		#Reaction rates
		"""All of the following rates are calculated by the equation
		r_ik = n_i*n_k/(rho*(1+delta_ik)*lambda_ik)
		where i,k define the elements and delta_ik =1 if i=k else 0.
		"""
		lambda_pp = rrs[0]/N_A*1E-6		#[m^3/s]
		lambda_pd = 1				#This one is actually unknown so it's just \
							#set to 1 since the reaction it is used in 
							#includes the number of deuterium atoms and \
							#that number is set to 0 earlier.
		lambda_33 = rrs[1]/N_A*1E-6		#[m^3/s]
		lambda_34 = rrs[2]/N_A*1E-6		#[m^3/s]
		lambda_7e = rrs[3]/N_A*1E-6		#[m^3/s]
		lambda_71prime  = rrs[4]/N_A*1E-6	#[m^3/s]
		#lambda_71  = Berylliym_Boron/N_A*1E-6	#[m^3/s]#Not used in this assigment
		#lambda_p14 = Nitrogen_Oxygen/N_A*1E-6	#[m^3/s]#Not used in this assigment

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
		#r_71  	 = n_Li*n_p	/rho*lambda_71	#Not used in this assigment
		#r_p14 	 = n_p*n_4	/rho*lambda_p14 #Not used in this assigment
		
		#This part makes sure that reaction rates which rely on other rates \
		#don't use elements that are not yet present.
		if (2./3)*r_33 + (1./3)*r_34 > r_pp:
			r_34 = (1./3)*r_pp
			r_33 = (2./3)*r_pp
		if r_7e > r_34:
			r_7e = r_34
		if r_71prime > r_7e:
			r_71prime = r_7e
		reactionrate_new = r_pp,r_pd,r_33,r_34,r_7e,r_71prime
		#print 'round:',i

		if array_equal(reactionrate_old,reactionrate_new) == True:
			counter +=1
		if counter == 20:
			counter = 0
			break
		reactionrate_old = reactionrate_new

		#Evolution of element abundances

		"""This part updates the number of Hydrogen atoms and Helium 4 atoms.
		We choose to neglect updating the rest of the elements since there \
		are so few of them, and thus their contribution is small."""
		d_npdt   = -n[2]**2.*lambda_pp 	- n[1]*n[2]*lambda_pd + n[3]**2.*\
			   lambda_33 - n[5]*n[2]*lambda_71prime
		d_nHe4dt = (n[3]**2.)/2.*lambda_33 - n[3]*n[4]*lambda_34 + 2.*n[5]\
			   *n[2]*lambda_71prime 
		n[2] 	+= dt*d_npdt 	    	#Number of Hydrogen atoms
		n[4] 	+= dt*d_nHe4dt   	#Number of Helium 4 atoms
	return reactionrate_new

#Solving the equation for the change in Luminosity

def Energy(T,rho):
	reactionrate = Reactionrates(T,rho)
	e_1  = reactionrate[0]*(Q_pp+Q_dp)
	e_2  = reactionrate[2]*Q_33
	e_3  = reactionrate[3]*Q_34
	e_4  = reactionrate[4]*Q_7e
	e_5  = reactionrate[5]*Q_71prime
	e    = e_1+e_2+e_3+e_4+e_5
	return e

#Defining functions for pressure,temperature and density
"""This part solves the equation for Pressure, assuming 
the equation of state to be that of an ideal gas"""

#Constants needed
mu 	= 1/(2.*X + 7./4*Y + 5./7*Z_74Be + 4./7*Z_73Li + 9./16*Z)
a 	= 4.*sigma/c

def Pressure(rho,T):
	#P_rad 	= a/.3*T**4		#Radiation pressure
	P_G 	= rho*k*T/(mu*m_u)    	#Gas pressure, assuming ideal gas
	P	= P_G #+ P_rad 		#Total Pressure
	return P

def Density(P,T):
	rho  = (P*mu*m_u)/(k*T)		#Ideal gas
	#rho = (m_u*mu)/(k*T)*(P - (a*T**4)/3.)	#Includes radiation pressure
	return rho

def Temperature(rho,P):
	T = (P*mu*m_u)/(k*rho)
	return T

"""
#Start of test area for Energy Production
#Test with core of the sun
rho     = 1.62E5    #Mass density of star [kg*m**-3]
T       = 1.57E7    #Temperature    [K]
X       = .7        #Hydrogen fraction
Y3      = 1E-10     #Helium 3 fraction
Y 	= .29 	    #Sum of Helium 3 and Helium 4 fractions
Y4      = Y-Y3      #Helium 4 fraction
Z_73Li  = 1E-7      #Lithium  fraction
Z_74Be  = 1E-7      #Berylium fraction
T   	= 1.57E7
rho 	= 1.62E5
reactionrate = Reactionrates(T,rho)
print "%.2e"%(reactionrate[0]*(Q_pp + Q_dp)*rho)
print "%.2e"%(reactionrate[2]*Q_33*rho)
print "%.2e"%(reactionrate[3]*Q_34*rho)
print "%.2e"%(reactionrate[4]*Q_7e*rho)
print "%.2e"%(reactionrate[5]*Q_71prime*rho)

dLdm = Energy(T_0,rho_0)
print "dL/dm = %.2e"%dLdm

P_0 = Pressure(rho_0,T_0)
print "P_0 = %.2e"%P_0
#End of test area for Energy Production
"""

#Function that reads the opacity table
lines = genfromtxt('opacity.txt')

def Opacity(T,rho):
	'''This function reads the opacity table'''
	#Convert R and T to the form used in opacity.txt
	R      = (rho*1E-3)/(T/1E6) 

	Rvalue = log10(R)
	Tvalue = log10(T)
	Rvalue = (round(2*Rvalue)/2.)
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
	kappa = (10.**(log10kappa/10.))	
	return kappa


#Solving the five equations
'''This part solves all the five equations'''

#Make empty lists
R = []
L = []
P = []
T = []
M = []
rho = []
epsilon = []

#Fill the arrays with the initial conditions
M.append(M_0)
R.append(R_0)
L.append(L_0)
T.append(T_0)
rho.append(rho_0)
#P.append(P_0)
P.append(Pressure(rho[0],T[0]))
#rho.append(Density(P[0],T[0]))
#T.append(Temperature(rho[0],P[0]))
epsilon.append(Energy(T[0],rho[0]))


#Simulate evolution by reducing mass.
i = 0     #Starts the counter
p = 0.01 #Max change allowed

while M[i] > 0 and R[i] > 0 and L[i] > 0 and T[i] > 0 and P[i] > 0:

    #Change in radius
    drdm 	= 1./(4*pi*R[i]**2*rho[i])
    dmr 	= p*R[i]/drdm    

    #change in Pressure
    dPdm 	= -(G*M[i])/(4*pi*R[i]**4)
    dmP 	= p*P[i]/dPdm        

    #Change in Luminosty
    dLdm 	= epsilon[i]		
    dmL 	= p*L[i]/dLdm    

    #Change in temperature
    kappa   	= Opacity(T[i],rho[i])		#Opacity
    dTdm    	= -3*kappa*L[i]/(256*pi**2*sigma*R[i]**4*T[i]**3)
    dmT 	= p*T[i]/dTdm

    #Change in density
    drhodm    	= (mu*m_u)/(k*T[i])*dPdm - (P[i]*mu*m_u)/(k*T[i]**2)*dTdm
    dmrho	= p*rho[i]/drhodm

    #Decide how small dm needs to be
    dmM =   10*p*M[i]
    dm  =   -min(abs(dmr),abs(dmP),abs(dmL),abs(dmT),abs(dmrho),abs(dmM))
    if abs(dm) < 1E3:
        dm = -1E3
    #Update values
    rho.append(rho[i] + drhodm*dm)
    M.append(M[i] + dm)			#Update mass
    R.append(R[i] + drdm*dm)		#Update Radius
    L.append(L[i] + dLdm*dm) 		#Update Luminosity
    P.append(P[i] + dPdm*dm)		#Update Pressure
    T.append(T[i] + dTdm*dm)		#Update Temperature
    epsilon.append(Energy(T[i],rho[i])) #Update energy output
    
    print 'M = %.3e, R = %.3e, L = %.3e, P = %.3e, T= %.3e, rho = %.3e, \
kappa = %.3e, epsilon = %.3e, dT = %.3e, dr = %.3e, dL = %.3e, dm = %.3e, \
dP = %.3e, drho = %.3e'%(M[i], R[i], L[i], P[i], T[i], rho[i], kappa, \
epsilon[i], dTdm*dm, drdm*dm, dLdm*dm, dm, dPdm*dm, drhodm*dm)


    i +=1   #increase counter

    #Print  information about which variable went to zero first.   
    if M[i] < 0:
        print 'Mass went negative'
	break
    if R[i] < 0:
        print 'Radius went negative'
	break
    if L[i] < 0:
        print 'Limuniosity went negative'
	break
    if T[i] < 0:
        print 'Temperature went negative'
	break
    if P[i] < 0:
        print 'Pressure went negative'
	break
  


#Scale values after initial conditions
"""for i in range(1,len(M)):
    M[i] = M[i]/M[0]
    
for i in range(1,len(R)):
    R[i] = R[i]/R[0]
    
for i in range(1,len(L)):
    L[i] = L[i]/L[0]
    
for i in range(1,len(P)):
    P[i] = P[i]/P[0]
    
for i in range(1,len(T)):
    T[i] = T[i]/T[0]
"""
#Plots of interesting things

figure()
plot(M,R)
title('Radius vs Mass')
xlabel('Mass[kg]')
ylabel('Radius[m]')

figure()
plot(M,P)
title('Pressure vs Mass')
xlabel('Mass[kg]')
ylabel('Pressure[Pa]')

figure()
plot(M,L)
title('Luminosity vs Mass')
xlabel('Mass[kg]')
ylabel('Luminosity[W]')

figure()
plot(M,T)
title('Temperature vs Mass')
xlabel('Mass[kg]')
ylabel('Temperature[K]')

figure()
plot(M,rho)
title('Density vs Mass')
xlabel('Mass[kg]')
ylabel('Density[kg/m^3]')

figure()
plot(M,epsilon)
title('Energyproduction vs Mass')
xlabel('Mass[kg]')
ylabel('Energy[J/kg]')

show()
