from scitools.all import *
#Values from text
E = array([5.,0.,0.])
m = 2.
q = 3.
dt = 10**-4
T  = 1.
n  = int(T/dt)
F = q*E

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
    a = F/m
    v[i+1] = v[i] + a*dt
    r[i+1] = r[i] + v[i+1]*dt

#Analytical solution
def R(t):
    return 15./4*t**2

#Plot movement in x-direction

plot(t,r[:,0],"r-") #Rx
legend("Analytical")
hold("on")

plot(t,R(t),"g-")
legend("Numerical")
xlabel("Time")
ylabel("Position")

hardcopy("opg1ab.png")


#We see that the numerical solution and the analytical lie
#perfectly on top of each other.
