from scitools.std import *

star = raw_input("name of datafile: ")

def pec_vel_func(wlobs):		# Calculates peculiar velocity with Doppler
    return c*( (wlobs*(10**-9)-wl) / wl )

t,wlobs,f = loadtxt(star,unpack=True) 	#Import time, observed wavelength and flux from file
c   = 2.99792458*10**8 			#Speed of light
wl = 656.3*10**(-9) 			#Reference wavelength

pec_vel = pec_vel_func(wlobs)		#Saves peculiar velocity 
pec_vel_mean = mean(pec_vel)		#Calculates mean peculiar velocity
rel_vel = pec_vel-pec_vel_mean		#Calculates relative velocity

#Plotting
figure(0)
plot(t,rel_vel)
title(star)
xlabel("Time")
ylabel("Relative velocity")
savefig(raw_input("Name of velocity figure: "))

figure(1)
plot(t,f)
xlabel("Time")
ylabel("Flux")
savefig(raw_input("Name of flux figure: "))
#wavelength results
#ja:  0,3,4,5,9
#nei: 1,2,7,8

#flux results
#ja:  0,3,4,
#nei: 1,2,5,6,7,8,9



#Measuring planet masses
G = 6.67384*10**-11	#Gravitational Constant

def mp(ms,vs,P):
	return (ms**(2./3)*vs*P**(1./3)) / ( (2*pi*G)**(1./3) )

ms3,ms0,ms4 = 0.5*(1.9891*10**30), 0.8*(1.9891*10**30), 1.8*(1.9891*10**30)
vs3,vs0,vs4 = 765, 225, 320
P3 ,P0 ,P4  = 170000,350000,1000000

mp3 = mp(ms3,vs3,P3)
mp0 = mp(ms0,vs0,P0)
mp4 = mp(ms4,vs4,P4)

print "Mass of planets"
print "Planet3: ",mp3
print "Planet0: ",mp0
print "Planet4: ",mp4

"""
thorae@rubin ~/privat/sem3/ast1100/oblig2 $ python Andreas_Ellewsen_Oblig2_del1.py 
Mass of planets
Planet3:  5.64127066095e+27
Planet0:  2.88747156539e+27
Planet4:  1.00057774805e+28
"""
