from scitools.std import *

def pec_vel_func(wlobs):
    return c*( (wlobs*(10**-9)-wl) / wl )

c  = 2.99792458*10**8 #m/s
wl = 656.3*10**(-9) #[m]

#Star 3
t,wlobs,f = loadtxt("star3.txt",unpack=True)	#Import info from file
pec_vel = pec_vel_func(wlobs)
pec_vel_mean = mean(pec_vel)
rel_vel = pec_vel-pec_vel_mean

t0=linspace(125000,131000,40)	#Estimated timeinterval around vrmax
vr=linspace(760,770,40)		#Estimated interval around vrmax
P = linspace(165000, 175000,40)	#Estimated interval with one period

best_t0    = t0[0]	
best_vr    = vr[0]
best_P     =  P[0]
best_delta = 10**60

#For loop for the function
for i in range(len(t0)):
	for j in range(len(vr)):
		for k in range(len(P)):
			rel_vel_model = vr[j]*cos(2*pi/P[k]*(t-t0[i]))
			delta = sum((rel_vel - rel_vel_model)**2)
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

#Star 0
#Exactly the same as Star3 is done
t,wlobs,f = loadtxt("star0.txt",unpack=True)	
pec_vel = pec_vel_func(wlobs)
pec_vel_mean = mean(pec_vel)
rel_vel = pec_vel-pec_vel_mean

t0=linspace(290000,304000,20)
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

#Star4
#Exactly the same as Star3 is done
t,wlobs,f = loadtxt("star4.txt",unpack=True)

pec_vel = pec_vel_func(wlobs)
pec_vel_mean = mean(pec_vel)
rel_vel = pec_vel-pec_vel_mean

t0=linspace(480000,540000,20)
vr=linspace(320,340,20)
P = linspace(1000000, 1200000,20)

best_t0    = t0[0]
best_vr    = vr[0]
best_P     =  P[0]
best_delta = 10**60

ape = 0
for i in range(len(t0)):
	for j in range(len(vr)):
		for k in range(len(P)):
			rel_vel_model = vr[j]*cos(2*pi/P[k]*(t-t0[i]))
			delta = sum((rel_vel - rel_vel_model)**2)
			if delta<best_delta:
				best_t0 = t0[i]
				best_vr = vr[j]
				best_P  =  P[k]
			ape+=1

figure(3)
plot(t,best_vr*cos(2*pi/best_P*(t-best_t0)))
hold("on")
plot(t,rel_vel)
xlabel("Time[s]")
ylabel("Velocity[m/s]")
title("star4")
