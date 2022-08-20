from scitools.std import *

#Values for timesteps
t_start = 0
t_stop  = 10.
dt = 0.001
n = int(t_stop/dt)

#Preparing arrays
t = linspace(t_start,t_stop,n+1)
x = zeros(len(t))
v = zeros(len(t))
a = zeros(len(t))

#Initial conditions
a[0] = 5
v[0] = 0
x[0] = 0

#Values in formula for acceleration
m= 80.
F= 400.
rho= 1.293
Cd= 1.2
A= 0.45
w= 0

#Euler method with plot
hold("on")
for k in range(n):
    v[k+1] = v[k] + dt*a[k]
    x[k+1] = x[k] + dt*v[k]
    a[k+1] = ( 2*F - rho*Cd*A*(v[k+1]-w)**2 ) / (2*m)
    if x[k+1]>100:
        print 'Reached %sm after %ss, with a velocity \
of %sm/s and an acceleration of %sm/s**2' %(x[k+1], t[k+1], v[k+1], a[k+1])
        plot(t,v, t,a, t,x, axis=[0,t[k+1],0,x[k+1]],legend=['v(t)','a(t)','x(t)'],xlabel= t, hardcopy='plotOppgaveE.png')
        break
