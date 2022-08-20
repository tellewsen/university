#Exercise E.45 p.727

#m*u'' + f(u') + s(u) = F(t)	, t>0
#u'' = ( F(t) - f(u') - s(u) ) / m

#u0 = u

#u[0]' = v[1]
#u[1]' = ( F(t) - f(u[1]) - s(u[0]) ) / m


def rhs(u,t):
    def f(dudt):	#friction
        return 0	#f(u')

    def s(u):	#spring
        return u	#s(u)

    def F(t):	#external
        return 0	#F(t)
    m = 1
    return [u[1], ( F(t) - f(u[1]) - s(u[0]) ) / m]

from scitools.std import *
from ODESolver import ODESolver,RungeKutta4,ForwardEuler
from RungeKutta2 import RungeKutta2
u0 = (1,0)
dt = pi/20
T = 2*pi
time_points = linspace(0,2*pi,T/float(dt))

plot(time_points, cos(time_points),'bo',
     time_points, -sin(time_points),'ro',
     legend= ('Exact( cos(t) )','Exact( -sin(t) )'))
hold('on')
for i in 0,1:
    for method in RungeKutta4,ForwardEuler,RungeKutta2:
        solver = method(rhs)
        solver.set_initial_condition(u0)
        u,t = solver.solve(time_points)
        figure(1)
        plot(t,u[:,i],legend= method.__name__)
        figure(2)
        plot(u[:,0],u[:,1], legend= method.__name__)


