#Exercise E.18 p.715

from scitools.std import *
from ODESolver import RungeKutta4, ForwardEuler

class Decay:				#Decay class
    def __init__(self,a,u0):		#Constructor
        self.a, self.u0 = float(a), u0

    def __call__(self,u,t):		#Decay function
        self.a = a
        return -a*u


u0 = 1					#Initial condition	given in exercise text
a  = log(2)/5600			#Given in exercisetext

carbon14time = Decay(a,u0)		#Carbon14 instance
dt = 500.				#Timesteps 		given in exercise text
T  = 20000				#Time in years 		given in exercise text
time_points = linspace(0,T,T/dt)	#Timestep array
Exact = exp(-a*T)			#Exact solution


print '%13s %10s %17s' %('Method', 'Final', 'Difference')	#A bad attempt at making this look good

for solver in ForwardEuler, RungeKutta4:		#Forloop for solving with the different methods
    problem = solver(carbon14time)			#Makes an instance with the method
    problem.set_initial_condition(carbon14time.u0)	#Sets the initial condition
    u, t = problem.solve(time_points)			#Solves with the method through time_points array
    plot(t,u, legend = solver.__name__)			#Plots that solution
    hold('on')	
    Final = u[-1]					#Saves last value in the solution array as Final
    Difference = Exact - Final				
    #Difference = abs(Exact - Final) Switch line above with this one if you dont like negative values
    print '%14s: %13.10f %13.10f' %(solver.__name__, Final, Difference) #Prints method,solution,diff from exact

"""
Terminal> python radioactive_decay.py 
       Method      Final        Difference
  ForwardEuler:  0.0774917253  0.0066270367
   RungeKutta4:  0.0841187917 -0.0000000297
"""
