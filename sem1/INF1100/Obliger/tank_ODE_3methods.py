
from scitools.std import *
from ODEsolver import *
class Problem:
    def __init__(self,h0,r,R,dt,T):
        self.h0,self.r,self.R,self.dt,self.T = h0,r,R,dt,T
        self.time_points = linspace(0, T, int(T/dt))

    def f(self,h,t):
        r,R = self.r, self.R
        g = 9.81
        if h<0:
            h = 0
        return -(float(r)/R)**2*sqrt(2*g*h)


def solve(problem,method):
    m = method(problem.f)
    m.set_initial_condition(problem.h0)
    h,t = m.solve(problem.time_points)
    plot(t,h, legend="%g  %s" %(problem.dt, method))
    

    
h0 = 1
r = 0.01
R = 0.2
T=200
problems = [Problem(h0,r,R,float(dt), T) for dt in sys.argv[1:]]
for problem in problems:
    for method in ForwardEuler,RungeKutta4,BackwardEuler:
        print "hei"
        solve(problem,method)
        hold('on')
r = raw_input(' ')
