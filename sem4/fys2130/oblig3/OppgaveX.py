# Oppgave X
from numpy import *
from matplotlib.pyplot import *

def diffEq(zi, vi, ti):
	ai = (-k*zi-b*vi)/m
	return ai

def rk(zstart, vstart, tstart):
	a1 = diffEq(zstart, vstart, tstart)
	v1 = vstart

	zhalf1 = zstart + v1*(dt/2.0)
	vhalf1 = vstart +a1*(dt/2.0)

	a2 = diffEq(zhalf1, vhalf1, tstart+dt/2.0)
	v2 = vhalf1

	zhalf2 = zstart + v2*(dt/2.0)
	vhalf2 = vstart +a2*(dt/2.0)

	a3 = diffEq(zhalf2, vhalf2, tstart+dt/2.0)
	v3 = v2

	zend = zstart + v3*dt
	vend = vstart + a3*dt

	a4 = diffEq(zend, vend, tstart+dt)
	v4 = vend

	amiddle = (1.0/6)*(a1+2*a2+2*a3+a4)
	vmiddle = (1.0/6)*(v1+2*v2+2*v3+v4)

	zend = zstart + vmiddle*dt
	vend = vstart + amiddle*dt
	
	return zend, vend
	
#Assignment starts here
#a)
#Underkritisk
#Constants
m = 0.1 #kg
k = 10.0  #N/m 
b = 0.1  #kg/s
T = 100.0
N = 10000

#Arrays
zlist = zeros(N)
vlist = zeros(N)
tlist = linspace(0,T,N)
dt = T/N

#Initial conditions
zlist[0] = 0.1 #m
vlist[0] = 0   #m/s
tlist[0] = 0   #s

#Rk4 integrator
for i in xrange(N-1):
	zlist[i+1], vlist[i+1]= rk(zlist[i], vlist[i], tlist[i])
plot(tlist,zlist)
xlabel('time[s]')
ylabel('vertical position[m]')
title('Damped oscillation (RK4)')
savefig("OppgaveXa.png")
figure()

#b)
#Underkritisk
#Constants
m = 100.0 #kg
b = 30.  #kg/s
k = 5*b**2/(4*m) #N/m 

#Arrays
zlist = zeros(N)
vlist = zeros(N)
tlist = linspace(0,T,N)
dt = T/N

#Initial conditions
zlist[0] = 0.1 #m
vlist[0] = 0   #m/s
tlist[0] = 0   #s

#Rk4 integrator
for i in range(N-1):
	zlist[i+1], vlist[i+1]= rk(zlist[i], vlist[i], tlist[i])
plot(tlist,zlist)
xlabel('time[s]')
ylabel('vertical position[m]')
title('Damped oscillation (RK4)')

#Kritisk

#Constants
m = 100.0 #kg
b = 30.  #kg/s
k = b**2/(4*m)  #N/m 

#Arrays
zlist = zeros(N)
vlist = zeros(N)

#Initial conditions
zlist[0] = 0.1 #m
vlist[0] = 0   #m/s
tlist[0] = 0   #s 

for i in range(N-1):
	zlist[i+1], vlist[i+1]= rk(zlist[i], vlist[i], tlist[i])
plot(tlist,zlist,'r-')

#Overkritisk
#Constants
m = 100.0 #kg
b = 30.  #kg/s
k = 0.5*b**2/(4*m)  #N/m 

#Arrays
zlist = zeros(N)
vlist = zeros(N)

#Initial conditions
zlist[0] = 0.1 #m
vlist[0] = 0   #m/s
tlist[0] = 0   #s 

for i in range(N-1):
	zlist[i+1], vlist[i+1]= rk(zlist[i], vlist[i], tlist[i])
plot(tlist,zlist,'c-')
legend(['Underkritisk','Kritisk','Overkritisk'])
savefig('OppgaveXb.png')
show()

