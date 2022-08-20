from scitools.std import *

#Constants
m1 = 6.4*10**23      	#[kg] mass of Mars
m2 = 1000.            	#[kg] mass of Satellite
G  = 6.67384*10**-11 	#[m^3kg^-1s^-2] Gravitational Constant

#timesteps
n  = 10**5
dt = 1. 				#[s]
time = n*dt			 	#[s]

#Arrays
r1 = zeros(n)  			#distance from origo to Mars
r2 = zeros(n)  			#distance from origo to satellite
r  = zeros(n)  			#distance between masses
rx1 = zeros(n) 			#position of Mars
ry1 = zeros(n) 			#position of Mars
rx2 = zeros(n) 			#position of Satellite
ry2 = zeros(n) 			#position of Satellite
ax1 = zeros(n)			#Acceleration Mars x
ay1 = zeros(n)			#Acceleration Mars y
ax2 = zeros(n)			#Acceleration Satellite x
ay2 = zeros(n)			#Acceleration Satellite y
vx1 = zeros(n)			#Velocity Mars x
vy1 = zeros(n)			#Velocity Mars y
vx2 = zeros(n)			#Velocity Satellite x
vy2 = zeros(n)			#Velocity Satellite y
F1x = zeros(n)
F1y = zeros(n)
F2x = zeros(n)
F2y = zeros(n)
#Initial conditions
ax1[0] = 0
ay1[0] = 0
vx1[0] = 0
vy1[0] = 0
ax2[0] = 0
ay2[0] = 0
vx2[0] = 0
vy2[0] = 1166 #[m/s]
rx1[0] = 0 					#[m^3] start postition of Mars 
ry1[0] = 0 					#[m^3] start postition of Mars
rx2[0] = float(10107+3400) 			#[m^3] start position of Satellite
ry2[0] = 0 					#[m^3] start postition of Satellite
r1[0] = sqrt(rx1[0]**2+ry1[0]**2)
r2[0] = sqrt(rx2[0]**2+ry2[0]**2)
r[0]  = sqrt((rx1[0]+rx2[0])**2+(ry1[0]+ry2[0]**2))

#Euler-Cromer for loop
for i in range(n-1):
	F1x[i+1] = (G*m1*m2)/(r[i]**3)*(rx1[i]+rx2[i]) #TROR FEILEN LIGGER HER
	F1y[i+1] = (G*m1*m2)/(r[i]**3)*(ry1[i]+ry2[i]) #OG HER
	F2x[i+1] = -F1x[i] #MULIG DETTE BLIR FOR ENKELT OG
	F2y[i+1] = -F1y[i] #OG DET HER DA SELVFLGELGI
	ax1[i+1] = F1x[i+1]/m1 #FIKS lengde som AU isteden forresten.. 
	ay1[i+1] = F1y[i+1]/m1 #virker som bedre benevning her
	ax2[i+1] = F2x[i+1]/m2
	ay2[i+1] = F2y[i+1]/m2
	vx1[i+1] = vx1[i] +ax1[i+1]*dt
	vy1[i+1] = vy1[i] +ay1[i+1]*dt
	vx2[i+1] = vx2[i] +ax2[i+1]*dt
	vy2[i+1] = vy2[i] +ay2[i+1]*dt
	rx1[i+1] = rx1[i] +vx1[i+1]*dt
	ry1[i+1] = ry1[i] +vy1[i+1]*dt
	rx2[i+1] = rx2[i] +vx2[i+1]*dt
	ry2[i+1] = ry2[i] +vy2[i+1]*dt
	r1[i+1]  = sqrt(rx1[i+1]**2+ry1[i+1]**2)
	r2[i+1]  = sqrt(rx2[i+1]**2+ry2[i+1]**2)
	r[i+1] = sqrt((rx1[i+1]+rx2[i+1])**2+(ry1[i+1]+ry2[i+1])**2)

plot(rx1,ry1,"-r")
axis([-10000000,10000000,-10000000,10000000])
xlabel("Meters[m]")
ylabel("Meters[m]")
hold("on")
plot(rx2,ry2,"-b")



#DETTE FUNKER IKKE SOM DET SKAL...................
