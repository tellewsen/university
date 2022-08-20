#Exercise E.19 p. 715

from scitools.std import *
from ODESolver import *
import sys
import numpy as np
class Problem:
    def __init__(self,h_0,r,R,dt,T):
        self.h_0 = h_0
        self.r  = r
        self.R  = R
        self.dt = dt
        self.T  = np.linspace(0,T,float(T)/dt)
    def __call__(self,h):
        if h < 0:
            h = 0
        return -(r/R)**2*np.sqrt(2*g*h)

h_0 = 1.
r   = 0.01
R   = 0.2
g   = 9.81
T   = 100.


# Given h0, r, R
problems = [Problem(h_0,r,R,float(dt),T) for dt in sys.argv[1:]]
for problem in problems:
    for method in ForwardEuler,BackwardEuler,RungeKutta4:
        time_points = np.linspace(0,T,T/float(dt))
        problem1.set_initial_condition = h_0
        problem1 = method.solve(time_points)

        hold('on')



"""
Jeg gir opp. Kan du sende meg en l0sning pA denne? Fatter ikke hvordan jeg skal fA det til A virke..
Blir bare mer og mer forvirra. Hjelper sikkert ikke pA saken at jeg ikke fikk til den forrige tank oppgaven heller.
"""
