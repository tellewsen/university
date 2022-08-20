from scitools.all import *

#Given values
m = 2
q = 3
B = [0,0,3]
r0 = [0,0,0]
v0 = [5,0,0]

#Timesteps
T = 5.
dt = 10**-4
n = int(T/dt)

#Arrays
v = zeros((n,3))
r = zeros((n,3))
a = zeros((n,3))
t = linspace(0,5,n)

#Initial conditions
v[0]= v0
r[0]= r0

#Euler-Cromer
for i in range(n-1):
    F = q*cross(v[i],B)    
    a[i+1] = F/m
    v[i+1] = v[i] + a[i+1]*dt
    r[i+1] = r[i] + v[i+1]*dt

#Plot position in each direction
figure(0)
plot(t,r[:,0],"r-")
legend("Rx")
hold("on")
plot(t,r[:,1],"g-")
legend("Ry")
plot(t,r[:,2],"b-")
legend("Rz")
hold("off")
hardcopy("opg2apos.png")

#Plot speed in each direction
figure(1)
plot(t,v[:,0],"r-")
legend("Vx")
hold("on")
plot(t,v[:,1],"g-")
legend("Vy")
plot(t,v[:,2],"b-")
legend("Vz")
hold("off")
hardcopy("opg2avel.png")

#Plot position in 3D
figure(2)
plot3(r[:,0],r[:,1],r[:,2])
hardcopy("opg2a3D.png")

#By looking at the plot, this period of the particle is approximately 1.4
print "Analytical solution: ",4*pi/9.
#Which is almost the same as the analytical, making this a pretty good model


