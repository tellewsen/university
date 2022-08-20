# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 21:59:33 2014

@author: Andreas Ellewsen
"""
from numpy import *
from matplotlib import * 
#This part is copied from wikipedia
def binomialCoefficient(n, k):
    if k < 0 or k > n:
        return 0
    if k > n - k: # take advantage of symmetry
        k = n - k
    if k == 0 or n <= 1:
    	return 1
    return binomialCoefficient(n-1, k) + binomialCoefficient(n-1, k-1)
#Part from Wikipedia ends here

#Defines functions
def Omega(N,q):
    return binomialCoefficient(N+q-1,q)

def P(Na,Nb,qa,qb):
    return Omega(Na,qa)*Omega(Nb,qb)/Omega(Na+Nb,q)
    
#Given values
Na = 2
Nb = Na
N = Na+Nb
q  = 6
qlist = linspace(0,q,q+1)
print(qlist)
#Calculate probabilities and print them in terminal
Prob = zeros(q+1)
for qa in qlist:
    Prob[qa] = P(Na,Nb,qa,q-qa)
    print('P(N=%s,qa=%s)= '%(N,qa),Prob[qa])

#Plot probability as function of energy units
plot(qlist,Prob)
xlabel('q[units of energy]')
ylabel('P(%s,q)'%N)
title('Probability P(qa) as function of qa for all qa')
show()