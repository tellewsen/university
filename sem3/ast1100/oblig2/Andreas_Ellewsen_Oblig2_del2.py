from scitools.std import *
#Constants
G = 6.67384*10**-11	#Gravitational Constant
c  = 2.99792458*10**8 	#m/s
wl = 656.3*10**(-9) 	#[m]


#Doppler function
def pec_vel_func(wlobs):
    return c*( (wlobs*(10**-9)-wl) / wl )



#Star 3
t,wlobs,f = loadtxt("star3.txt",unpack=True)	#Import info from file
pec_vel = pec_vel_func(wlobs)			#Peculiar velocity from Doppler
pec_vel_mean = mean(pec_vel)			#Mean velocity of Peculiar velocity
rel_vel = pec_vel-pec_vel_mean			#Relative velocity

t0=linspace(125000,131000,40)			#Estimated timeinterval around vrmax
vr=linspace(760,770,40)				#Estimated interval around vrmax
P = linspace(165000, 175000,40)			#Estimated interval with one period

best_t0    = t0[0]	
best_vr    = vr[0]
best_P     =  P[0]
best_delta = 10**60

#For-loop for the functions
for i in range(len(t0)):
	for j in range(len(vr)):
		for k in range(len(P)):
			rel_vel_model = vr[j]*cos(2*pi/P[k]*(t-t0[i]))	#Relative velocity model
			delta = sum((rel_vel - rel_vel_model)**2)	#Delta function
			if delta<best_delta:
				best_t0 = t0[i]
				best_vr = vr[j]
				best_P  =  P[k]

#Plotting
figure(1)
plot(t,best_vr*cos(2*pi/best_P*(t-best_t0)))	#t vs model with best variables
hold("on")
plot(t,rel_vel)					#t vs observed velocity
xlabel("Time[s]")
ylabel("Velocity[m/s]")
title("Star3")
savefig("Star3.png")

#Measuring planetmass3
vs3 = best_vr
P3  = best_P

#Star 0
#Exactly the same as Star3 is done
t,wlobs,f = loadtxt("star0.txt",unpack=True)	
pec_vel = pec_vel_func(wlobs)
pec_vel_mean = mean(pec_vel)
rel_vel = pec_vel-pec_vel_mean

t0=linspace(290000,295000,20)
vr=linspace(205,225,20)
P = linspace(345000, 350000,20)

best_t0    = t0[0]
best_vr    = vr[0]
best_P     =  P[0]
best_delta = 10**60

for i in range(len(t0)):
	for j in range(len(vr)):
		for k in range(len(P)):
			rel_vel_model = vr[j]*cos(2*pi/P[k]*(t-t0[i]))
			delta = sum((rel_vel - rel_vel_model)**2)
			if delta<best_delta:
				best_t0 = t0[i]
				best_vr = vr[j]
				best_P  =  P[k]

figure(2)
plot(t,best_vr*cos(2*pi/best_P*(t-best_t0)))
hold("on")
plot(t,rel_vel)
xlabel("Time[s]")
ylabel("Velocity[m/s]")
title("Star0")
savefig("Star0.png")

#Measuring planetmass0
vs0 = best_vr
P0  = best_P

#Star4
#Exactly the same as Star3 is done
t,wlobs,f = loadtxt("star4.txt",unpack=True)

pec_vel = pec_vel_func(wlobs)
pec_vel_mean = mean(pec_vel)
rel_vel = pec_vel-pec_vel_mean

t0=linspace(460000,580000,20)
vr=linspace(300,380,20)
P = linspace(850000, 1200000,20)

best_t0    = t0[0]
best_vr    = vr[0]
best_P     =  P[0]
best_delta = 10**60

for i in range(len(t0)):
	for j in range(len(vr)):
		for k in range(len(P)):
			rel_vel_model = vr[j]*cos(2*pi/P[k]*(t-t0[i]))
			delta = sum((rel_vel - rel_vel_model)**2)
			if delta<best_delta:
				best_t0 = t0[i]
				best_vr = vr[j]
				best_P  =  P[k]
				best_delta = delta
figure(3)
plot(t,best_vr*cos(2*pi/best_P*(t-best_t0)))
hold("on")
plot(t,rel_vel)
xlabel("Time[s]")
ylabel("Velocity[m/s]")
title("star4")
savefig("Star4.png")

#Measuring planetmass4

vs4 = best_vr
P4  = best_P

#Planet mass measuring function
def mp(ms,vs,P):
	return (ms**(2./3)*vs*P**(1./3)) / ( (2*pi*G)**(1./3) )

#Mass of stars
ms3,ms0,ms4 = 0.5*(1.9891*10**30), 0.8*(1.9891*10**30), 1.8*(1.9891*10**30)
mp3,mp0,mp4 = mp(ms3,vs3,P3), mp(ms0,vs0,P0), mp(ms4,vs4,P4)

print "Mass of planets"
print "Planet3: ",mp3
print "Planet0: ",mp0
print "Planet4: ",mp4

"""
thorae@rubin ~/privat/sem3/ast1100/oblig2 $ python Andreas_Ellewsen_Oblig2_del2.py 
Mass of planets
Planet3:  5.73327274832e+27
Planet0:  2.88747156539e+27
Planet4:  1.06311385731e+28
"""
