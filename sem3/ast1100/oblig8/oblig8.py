from scitools.all import *
import random

#Constants
m = 1.67E-27		#Mass of hydrogen atom (one proton mass)
T = 6000.		#Temperature of the sun
k = 1.38E-23		#Boltzmanns constant
n = 10**5		#Iterations
sigma = sqrt(k*T/m)	#Standard deviation

#Velocities and Kinetic Energy for each particle
K = zeros(n)
V = zeros((n,3))
for i in range(n):
    V[i,0] = random.gauss(0,sigma)
    V[i,1] = random.gauss(0,sigma)
    V[i,2] = random.gauss(0,sigma)
    K[i]   = 1./2*m*(V[i,0]**2+V[i,1]**2+V[i,2]**2)

#Calculate mean kinetic energy for a particle both ways
numK  = mean(K)			#Numerical
anaK  = 3./2*k*T		#Analytical
Error = sqrt((numK-anaK)**2)	#Difference between the two

#Print results
print "Numerical  : ",numK
print "Analytical : ",anaK
print "Error      : ",Error
print "Error in %% :  %s%%" %(Error/anaK*100)

"""
thorae@safir ~/privat/sem3/ast1100/oblig8 $ python oblig8.py 
scitools.easyviz backend is gnuplot
Numerical  :  1.24115652859e-19
Analytical :  1.242e-19
Error      :  8.43471411701e-23
Error in % :  0.0679123519888%
"""
