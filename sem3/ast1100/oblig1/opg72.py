from scitools.std import *

#Constants
G  = 6.67384*10**-11		#[m^3 kg^-1 s^-2]  Gravitational Constant
m1 = 6.4*10**23			#[kg] 		Mass of Mars
m2 = 100.			#[kg] 		Mass of Beagle2
k  = 0.00016 			#[kg/s]		friction constant
dt = 1.				#[s^1]		timestep

#Position lists and initial position
x0_2=-3698000.			#[m]
y0_2= 0.			#[m]
r   = sqrt((x0_2)**2+(y0_2)**2)	#[m]

vx2 = 0				#[m/s]
vy2 =-4000. 			#[m/s]
ax2 = 0.			#[m/s**2]
ay2 = 0.			#[m/s**2]

x2  = []			#[m]
y2  = []			#[m]
x2.append(x0_2)
y2.append(y0_2)

#Counter
i = 0 
#Euler-Cromer while-loop
while r>3400000:
	#acceleration
	ax2 = -G*m1/r**3*x2[i] - k*vx2/m2
	ay2 = -G*m1/r**3*y2[i] - k*vy2/m2
	#velocity
	vx2 = vx2 + ax2*dt
	vy2 = vy2 + ay2*dt
	#position
	x2.append(x2[i] + vx2*dt)
	y2.append(y2[i] + vy2*dt)
	r =sqrt(x2[i+1]**2+y2[i+1]**2)
	i+=1

#Plotting
array(x2)
array(y2)
Mars=linspace(0,2*pi,1001)

plot(x2,y2,"b-")
axis("equal")
xlabel("Meters[m]")
ylabel("Meters[m]")
hold("on")

plot(3400000*cos(Mars),3400000*sin(Mars),"r--")
legend("Position of Beagle2","Mars")
title("Landing of Beagle2")
savefig("opg72.png")

#Can conclude that Beagle2 was supposed to study the rocks close to equator
