from scitools.all import *

#Given values
m = 2
q = 3
B = [0,0,3]
r0 = [0,0,0]
v0 = [5,0,2]

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

#Plots
plot3(r[:,0],r[:,1],r[:,2])
legend("Position")
hardcopy("opg2d3D.png")
