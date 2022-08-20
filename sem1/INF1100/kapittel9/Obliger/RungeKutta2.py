#Exercise E.22 p.714
from scitools.std import *
from ODESolver import *
class RungeKutta2(ODESolver):
    def advance(self):
        u, f, k, t = self.u, self.f, self.k, self.t
        dt = t[k+1] - t[k]
        K1 = dt*f(u[k],t[k])
        K2 = dt*f(u[k] +0.5*K1, t[k] + 0.5*dt)
        u_new = u[k] + K2
        return u_new


if __name__ == "__main__":
    # Test Problem
    def f(u,t):
        return u

    u0= 1
    time_points = linspace(0,2*pi,100)
    solver = RungeKutta2(f)
    solver.set_initial_condition(u0)
    u,t = solver.solve(time_points)

    Exact = exp(t)
    Difference = abs(u-Exact)

    #Plot of Test Problem
    plot(t,Difference,
        title="Difference RungeKutta2 vs Exact",
        xlabel= "t", 
        ylabel="Difference")

"""
Terminal> python RungeKutta2.py 
"""
