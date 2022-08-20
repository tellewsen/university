"""
Term Project 0
The first term project involves extracting the energy production at the
center of a star, given a temperature and density. The project involves
information contained in Chaps.1-3.

Created by Andreas Ellewsen
Started 18.02.15
"""

from math import exp

#Constants
m_u = 1.66053892173E-27 	#[kg]	Atomic mass unit
N_A = 6.0221413E23      	#Avogadros number
MeVtoJ = 1.60217657E-13 	#[J] Conversion factor from MeV to Joules

#Initial conditions
rho     = 1.62E5    #Mass density of star [kg*m**-3]
T       = 1.57E7    #Temperature    [K]
X       = .7        #Hydrogen fraction
Y3      = 1E-10     #Helium 3 fraction
Y 	   = .29	    #Sum of Helium 3 and Helium 4 fractions
Y4      = Y-Y3      #Helium 4 fraction
Z_73Li  = 1E-7      #Lithium  fraction
Z_74Be  = 1E-7      #Berylium fraction

#Q values

Q_pp = (.15+1.02)*MeVtoJ	#[J]
Q_dp = (5.49)*MeVtoJ		#[J]
Q_33 = (12.86)*MeVtoJ		#[J]
Q_34 = (1.59)*MeVtoJ		#[J]
Q_7e = (.05)*MeVtoJ		#[J]
Q_71prime = (17.35)*MeVtoJ	#[J]

#Atom numbers
n_e = rho/(2.*m_u)*(1+X)    #Number of electrons
n_d = 0			    #Number of deuterium atoms (assumed creation=destruction => 0)
n_p = rho*X/m_u     	    #Number of Hydrogen atoms
n_3 = rho*Y3/(3.*m_u)       #Number of Helium 3 atoms
n_4 = rho*Y4/(4.*m_u)       #Number of Helium 4 atoms
n_Li  = rho*Z_73Li/(7.*m_u) #Numver of Beryllium 7 atoms
n_Be  = rho*Z_74Be/(7.*m_u) #Number of Lithium 7 atoms


#Reactions
"""In this part the reaction rates of the processes in the PPI and PPII chain are calculated. We neglect the other processes since we choose to look at a star with low temperature at it's core, and in that case those processes are very small compared to PPI and PPII."""

T9 = T/1E9
T9_star1 = T9/(1+4.95E-2*T9)
T9_star2 = T9/(1+.759*T9)

Lithium_proton = (1.096E9*T9**(-2./3)*exp(-8.472*T9**(-1./3)) - 4.830E8*T9_star2**(5./6)*T9**(-3./2)*exp(-8.472*T9_star2**(-1./3)) + 1.06E10*T9**(-3./2)*exp(-30.442*T9**(-1.)))
Hydrogen_Deuterium  = 4.01E-15*T9**(-2./3)*exp(-3.380*T9**(-1./3))*(1 +
                      .123*T9**(1./3) + 1.09*T9**(2./3) + .938*T9)
Helium3_proton       = 6.04E10*T9**(-2./3)*exp(-12.276*T9**(-1./3))*(1 + 
                      .034*T9**(1./3) - .522*T9**(2./3)-.124*T9 + 
                      .353*T9**(4./3)+ .213*T9**(-5./3))
Helium3_Beryllium   = 5.61E6*T9_star1**(5./6)*T9**(-3./2)*exp(-12.826*T9_star1**(-1./3))
Beryllium_Lithium   = 1.34E-10*T9**(-1./2)*(1 - .537*T9**(1./3)  + 3.86*T9**(2./3) +  .0027*T9**-1*exp(2.515E-3*T9**-1))
Lithium_proton	    = 	1.096E9*T9**(-2./3)*exp(-8.472*T9**(-1./3)) - 4.830E8*T9_star2**(5./6)*T9**(-3./2)*exp(-8.472*T9_star2**(-1./3)) + 1.06E10*T9**(-3./2)*exp(-30.442*T9**-1)
"""
"These reactions are not used in this assignment, but may be useful at a later point, so I've kept them here."
Berylliym_Boron	    = 	3.11E5*T9**(-2./3)*exp(-10.262*T9**(-1./3))+2.53E3*T9**(-3./2)*exp(-7.306*T9**-1)
Nitrogen_Oxygen	    = 	4.9E7*T9**(-2./3)*exp(-15.228*T9**(-1./3)-
			.092*T9**2)*(1+.027*T9**(1./3)-.778*T9**(2./3)-.149*T9+.261*T9**(4./3)+
			.127*T9**(5./3)) +2.37E3*T9**(2./3)*exp(-3.011*T9**-1) + 2.19E4*exp(-12.53*T9**-1)
"""
#Reaction rates
"""All of the following rates are calculated by the equation
rik = n_i*n_k/(rho*(1+delta_ik)*lambda_ik)
where i,k define the elements and delta_ik =1 if i=k else 0.
"""

lambda_pp = Hydrogen_Deuterium/N_A*1E-6		#[m^3/s]
r_pp 	  = n_p**2/(rho*2)*lambda_pp		

lambda_pd  = 1					#This one is actually unknown so it's just set to 1 since the reaction it is used in 
						#includes the number of deuterium atoms and that number is set to 0 earlier.

r_pd  = r_pp					#Assume this reaction happens instantly such that it happens at the same rate as the elements it needs become available.
						#Thus it must be the same as the reaction that creates those elements (in this case Deuterium)
 
lambda_33 = Helium3_proton/N_A*1E-6		#[m^3/s]
r_33 	  = n_3**2/(rho*2)*lambda_33		

lambda_34 = Helium3_Beryllium/N_A*1E-6		#[m^3/s]
r_34 	  = n_3*n_4/rho*lambda_34

lambda_7e = Beryllium_Lithium/N_A*1E-6		#[m^3/s]
r_7e 	  = n_Be*n_Li/rho*lambda_7e

lambda_71prime  = Lithium_proton/N_A*1E-6	#[m^3/s]
r_71prime       = n_Li*n_p/rho*lambda_71prime

"""
"These reactions are not used in this assignment, but may be useful at a later point, so I've kept them here."

lambda_71  = Berylliym_Boron/N_A*1E-6		#[m^3/s]
r_71 	   = n_Li*n_p/rho*lambda_71

lambda_p14 = Nitrogen_Oxygen/N_A*1E-6		#[m^3/s]
r_p14 	   = n_p*n_4/rho*lambda_p14
"""

#Evolution of element abundances
"""In this segment I make a simplification by assuming that there
 is no evolution in the density of any atomic species except Hydrogen
 and Helium 4."""

#Time evolution
"""This segement takes care of the time evolution of the reactions. By doing this we hope to reach some kind of equilibrium over time, such that all the reaction rates stabilize at some level. This is achieved remarkably fast, as can be seen in the choice of time T and steps N."""

Time  = 1		#Total time to simulate [s]
N  = 10**4	#Number of steps
dt = Time/N	#Timesteps [s]

for i in range(0,N-1):
	r_pp 	 = n_p**2	/(2*rho)*lambda_pp		
	r_33 	 = n_3**2	/(2*rho)*lambda_33		
	r_34 	 = n_3*n_4	/rho*lambda_34
	r_71prime= n_Li*n_p	/rho*lambda_71prime
	#r_71  	 = n_Li*n_p	/rho*lambda_71	#Not used in this assigment as mentioned above
	#r_p14 	 = n_p*n_4	/rho*lambda_p14 #Not used in this assigment as mentioned above
	
	#This part makes sure that reaction rates which rely on other rates don't use elements that are not yet present.
	if (2./3)*r_33 + (1./3)*r_34 > r_pp:
		r_34 = (1./3)*r_pp
		r_33 = (2./3)*r_pp
	if r_7e > r_34:
		r_7e = r_34
	if r_71prime > r_7e:
		r_71prime = r_7e
	
	#This part updates the number of Hydrogen atoms and Helium 4 atoms.
	#We choose to neglect updating the rest of the elements since there are so few of them, and thus
	#their contribution is small.
	d_npdt   = -n_p**2*lambda_pp 	- n_d*n_p*lambda_pd + n_3**2*lambda_33 - n_Li*n_p*lambda_71prime
	d_nHe4dt =  (n_3**2.)/2.*lambda_33 - n_3*n_4*lambda_34 + 2.*n_Li*n_p*lambda_71prime 
	n_p 	+=  dt*d_npdt 	    	#Number of Hydrogen atoms
	n_4 	+=  dt*d_nHe4dt   	#Number of Helium 4 atoms


#Testing against sanity test
"""Should get the following values:
all in units Jm-3s-1
r_pp(Q_pp + Q_dp)rho 	= 4.04E2
r_33Q_33rho 		= 1.08E-6 
r_34Q_34rho		= 3.66E-5	
r_7eQ_7erho		= 1.15E-6
r_71Q_71rho		= 4.00E-4
"""


test1 = r_pp	 *(Q_pp + Q_dp)*rho
test2 = r_33	 *Q_33*rho
test3 = r_34	 *Q_34*rho
test4 = r_7e	 *Q_7e*rho
test5 = r_71prime*Q_71prime*rho
print "Test:"
print "%.2e"%(test1)
print "%.2e"%(test2)
print "%.2e"%(test3)
print "%.2e"%(test4)
print "%.2e"%(test5)
