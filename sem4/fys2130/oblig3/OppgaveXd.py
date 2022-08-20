from numpy import *
from matplotlib.pyplot import *

def diffEq(z, v, t,omegaf,omega0):
	a = -omega0**2*z - (b/m)*v + (F/m)*cos(omegaf*t)
	return a

def rk(zstart, vstart, tstart,omegaf,omega0):
	a1 = diffEq(zstart, vstart, tstart,omegaf,omega0)
	v1 = vstart

	zhalf1 = zstart + v1*(dt/2.0)
	vhalf1 = vstart +a1*(dt/2.0)

	a2 = diffEq(zhalf1, vhalf1, tstart+dt/2.0,omegaf,omega0)
	v2 = vhalf1

	zhalf2 = zstart + v2*(dt/2.0)
	vhalf2 = vstart +a2*(dt/2.0)

	a3 = diffEq(zhalf2, vhalf2, tstart+dt/2.0,omegaf,omega0)
	v3 = v2

	zend = zstart + v3*dt
	vend = vstart + a3*dt

	a4 = diffEq(zend, vend, tstart+dt,omegaf,omega0)
	v4 = vend

	amiddle = (1.0/6)*(a1+2*a2+2*a3+a4)
	vmiddle = (1.0/6)*(v1+2*v2+2*v3+v4)

	zend = zstart + vmiddle*dt
	vend = vstart + amiddle*dt
	
	return zend, vend
	
#Assignment starts here

#Constants
m = 0.10	#kg
k = 10.0	#N/m 
b = 0.04	#kg/s
F = 0.1		#N
T = 100.0
N = 1000

omega0 = sqrt(k/m)

f0 = omega0/(2*pi)
#Arrays
zlist = zeros(N)
vlist = zeros(N)
tlist = linspace(0,T,N)
omegaF = sqrt(omega0**2 - (b**2)/(2*m**2))
Elist = zeros(N)
flist = zeros(N)
dt = T/N

def E(z,v):
	E = 0.5*k*z**2 + 0.5*m*v**2
	return E

#Initial conditions
zlist[0] = 0 	#m
vlist[0] = 0    #m/s
tlist[0] = 0    #s


#Rk4 integrator
for j in range(0,400):
	omega = omegaF*j*0.005
	for i in xrange(N-1):
		zlist[i+1], vlist[i+1]= rk(zlist[i], vlist[i], tlist[i], omega,omega0)
	Elist[j+1] = mean(E(zlist[0.75*N:N],vlist[0.75*N:N]))
	flist[j+1] = omega/(2*pi)
plot(flist,Elist)
xlabel('Frequency[Hz]')
ylabel('Energy[J]')
title('Time averaged energy (RK4)')
grid('on')
savefig('OppgaveXd.png')
show()

