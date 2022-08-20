#Exercise E.21 p.714

from scitools.std import *

def RungeKutta2(f,u0,T,n):
    """Solve u'=f(u,t), u(0) = u0, with n steps until t= T."""
    t = zeros(n+1)
    u = zeros(n+1) #u[k] is the solution at time t[k]
    u[0] = u0
    t[0] = 0
    dt = T/float(n)    
    for k in range(n):		#RungeKutta2 advance
        K1 = dt*f(u[k],t[k])
        K2 = dt*f(u[k] +0.5*K1, t[k] + 0.5*dt)
        t[k+1] = t[k] + dt
        u[k+1] = u[k] + K2
    return u,t

# Test Problem
def f(u,t):
    return u

u, t = RungeKutta2(f, u0=1, T=2*pi, n=100)
Exact = exp(t)
Difference = abs(u-Exact)

#Plot of Test Problem
plot(t,Difference,
    title="Difference RungeKutta2 vs Exact",
    xlabel= "t", 
    ylabel="Difference")


"""
Terminal> python RungeKutta2_func.py 
"""
