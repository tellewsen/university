from scitools.all import *
from numpy.linalg import norm
#Constants
G  = 6.67E-11			#[m^3 kg^-1 s^-2]  	Gravitational Constant
m1 = 5.97219E24			#[kg] 			Mass of Earth
m2 = 1.9891E30			#[kg] 			Mass of The Sun
m3 = 7.34767309E22		#[kg]			Mass of The Moon

AU = 1.496E11 			#[m/AU]			1AU in meters
dt = 400.			#[s]			timestep
T  = 60*60*24*365
n  = int(T/dt)			#[]			iterations
print "Iterations: ",n

#Useful stuff
REarth     = 6371000.		#[m]			Radius Earth
RMoon      = 1737100.		#[m]			Radius Moon
DEarthMoon = (406700000+356400000)/2.		#[m]			Disctance Earth->Moon
OrbitEarth = 29800		#[m/s]			Orbitspeed Earth->Sun
OrbitMoon  = 1022.		#[m/s]			Orbitspeed Moon->Earth

#Arrays
a1 = zeros((n,3))
a2 = zeros((n,3))
a3 = zeros((n,3))

v1 = zeros((n,3))
v2 = zeros((n,3))
v3 = zeros((n,3))

r1 = zeros((n,3))
r2 = zeros((n,3))
r3 = zeros((n,3))

#Initial Conditions

v1[0] = [ 0,OrbitEarth,0]				#[m/s]
v2[0] = [ 0, 0 ,0]					#[m/s]
v3[0] = [ OrbitMoon,OrbitEarth,0]			#[m/s]

r1[0] = [ AU,   0,0]					#[m]
r2[0] = [ 0,    0,0]					#[m]
r3[0] = [ AU, REarth + DEarthMoon + RMoon ,0]	#[m]

#Gravity function
def gravity(mass1,mass2,r):
	return -G*mass1*mass2/(norm(r)**3) * r

#Euler-Cromer loop
counter= 0
percent= 0
for i in range(n-1):
	r12  	= r2[i]-r1[i]	#distance between planet and small star
	r13  	= r3[i]-r1[i]	#distance between planet and large star
	r23 	= r3[i]-r2[i]	#distance between small star and large star

	F12 	= gravity(m1,m2,r12)
	F13	= gravity(m1,m3,r13)
	F23	= gravity(m2,m3,r23)
	F21 	= -F12
	F31	= -F13
	F32 	= -F23

	a1[i+1] =  (F31 + F21)/m1
	a2[i+1] =  (F12 + F32)/m2
	a3[i+1] =  (F13 + F23)/m3

	v1[i+1] = v1[i] + a1[i+1]*dt
	v2[i+1] = v2[i] + a2[i+1]*dt
	v3[i+1] = v3[i] + a3[i+1]*dt

	r1[i+1] = r1[i] + v1[i+1]*dt
	r2[i+1] = r2[i] + v2[i+1]*dt
	r3[i+1] = r3[i] + v3[i+1]*dt

	counter +=1
	if counter==n/10.:
		percent +=10
		print "%.0f%% done"%percent
		counter = 0

r1AU = r1/AU
r2AU = r2/AU
r3AU = r3/AU
#Plotting 2D
plot(r1AU[:,0],r1AU[:,1],"r-")
legend("Earth")
hold("on")
plot(r2AU[:,0],r2AU[:,1],"g-")
legend("Sun")
plot(r3AU[:,0],r3AU[:,1],"b-")
legend("Moon")
xlabel("X-axis[AU]")
ylabel("Y-axis[AU]")
raw_input('Press enter to end program')
