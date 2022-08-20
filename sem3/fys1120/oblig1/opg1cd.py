from scitools.all import *
#Values from text
E = array([1.,2.,-5.])
m = 2.
q = 3.

F = q*E
dt = 10**-4
T  = 1.
n  = int(T/dt)

#Arrays
t = linspace(0,1,n)
a = zeros((n,3))
v = zeros((n,3))
r = zeros((n,3))
#Initial conditions
r0 = [0,0,0]
v0 = [0,0,0]

#Euler-Cromer
for i in range(n-1):
    a[i+1] = F/m
    v[i+1] = v[i] + a[i+1]*dt
    r[i+1] = r[i] + v[i+1]*dt

#Plot position in each direction

figure(0)
plot(t,r[:,0],"r-") #Rx
n*0.4,legend("Rx")
hold("on")
plot(t,r[:,1],"g-") #Ry
legend("Ry")
plot(t,r[:,2],"b-") #Rz
legend("Rz")
xlabel("Time")
ylabel("Position")
hold("off")
hardcopy("opg1cpos.png")

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
hardcopy("opg1cvel.png")


x = []
y = []
z = []
for i in n*0.1,n*0.2,n*0.3,n*0.4,n*0.5,n*0.6,n*0.7,n*0.8,n*0.9,n-1:
    x.append(r[i,0])
    y.append(r[i,1])
    z.append(r[i,2])

#Plot position in 3D
figure(2)
plot3(r[:,0],r[:,1],r[:,2],'c-')
legend("Position")
hold("on")
plot3(x,y,z,'ro')
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
legend("Position at times with same dt between them")

hardcopy("opg1d.png")
#We see that the numerical solution and the analytical lie
#perfectly on top of each other.
