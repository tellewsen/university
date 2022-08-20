from scitools.std import *

#Values for timesteps
t_start = 0
t_stop  = 20.
dt      = 0.001
n       = int(t_stop/dt)

#Preparing arrays
t  = linspace(t_start,t_stop,n+1)
x  = zeros(len(t))
v  = zeros(len(t))
a  = zeros(len(t))
Fc = zeros(len(t)) 
Fv = zeros(len(t))
A  = zeros(len(t))
D  = zeros(len(t))
F  = zeros(len(t))

#Initial conditions
a[0]  = 5.
#v[0]  = 0.
#x[0]  = 0.
#A[0]  = 0.
#Fc[0] = 0.
#Fv[0] = 0.
#D[0]  = 0.
#F[0]  = 0.


#Given values
m  = 80. #Input Usain Bolt's mass for his 2008 world record (94kg)
rho= 1.293
Cd = 1.2
A  = 0.45
w  = 0.
fc = 488.
tc = 0.67
fv = 25.8

#Euler method with plot
hold("on")
for k in range(n):
    v[k+1]  = v[k] + dt*a[k]
    x[k+1]  = x[k] + dt*v[k]
    D[k+1]  = (1./2)*A*(1 - 0.25*exp(-(t[k+1]/(tc)**2)))*rho*Cd*(v[k+1]-w)**2
    Fc[k+1] = fc*exp(-((t[k+1])/tc)**2)
    Fv[k+1] = fv*v[k+1]    
    F[k+1]  = 400.
    a[k+1]  = (F[k+1]+Fc[k+1]-Fv[k+1]-D[k+1])/m
    if x[k+1]>100:
	print 'Reached %sm after %ss, with a velocity of %sm/s and an acceleration of %sm/s**2' %(x[k+1], t[k+1], v[k+1], a[k+1])
	figure(0)
	plot(t,x, t,v, t,a, axis=[0,t[k],0,x[k]], legend=['x(t)','v(t)','a(t)'], xlabel='t', ylabel='x[m],v[m/s],a[m/s**2]', hardcopy="plotOppgaveI.png")
	figure(1)
	plot(t,F, t,Fc, t,Fv, t,D, axis=[0,t[k+1],-10,10+Fc[1]],legend=['F(t)','Fc(t)','Fv(t)','D(t)'],xlabel='t', ylabel='Force', hardcopy="plotOppgaveK.png")
        break
