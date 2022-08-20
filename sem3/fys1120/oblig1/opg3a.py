from scitools.all import *
#Values from text
B = array([0,0,1])
m = 1.
q = 1.
dt = 10**-4
T  = 50.
n  = int(T/dt)
rD = 2.6
omega = q/m * sqrt(B[0]**2+B[1]**2+B[2]**2)
#Arrays
t = linspace(0,50,n)
v = zeros((n,3))
r = zeros((n,3))

#Initial conditions
r0 = [0,0,0]
v0 = [0,0,0]

#Euler-Cromer

for i in range(n-1):
    if -0.1 <= r[i,0] <=0.1:
        E = array([cos(omega*t[i]),0,0])
    else:
        E = 0
    F = q*E + q*cross(v[i],B)
    a = F/m
    v[i+1] = v[i] + a*dt
    r[i+1] = r[i] + v[i+1]*dt

#Plot position in each direction
figure(0)
plot(t,r[:,0],"r-") #Rx
legend("Rx")
hold("on")
plot(t,r[:,1],"g-") #Ry
legend("Ry")
plot(t,r[:,2],"b-") #Rz
legend("Rz")
xlabel("Time")
ylabel("Position")
hold("off")
hardcopy("opg3apos.png")

#Plot speed in each direction
figure(1)
plot(t,v[:,0],"r-") #Vx
legend("Vx")
hold("on")
plot(t,v[:,1],"g-") #Vy
legend("Vy")
plot(t,v[:,2],"b-") #Vz
legend("Vz")
xlabel("Time")
ylabel("Speed")
hold("off")
hardcopy("opg3avel.png")

#Plot position in 3D
figure(2)
plot(r[:,0],r[:,1])
legend("Position")
hardcopy("opg3a3D.png")
print "Just wait for it to finish, this takes some time..."

