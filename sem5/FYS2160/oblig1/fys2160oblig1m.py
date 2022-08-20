# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 22:45:18 2014

@author: Andreas Ellewsen
"""
from numpy import *
from matplotlib import * 
import operator
import math

#Defines functions
def Omega(N,q):
    return math.factorial(N+q-1)/math.factorial(q)/math.factorial(N-1)

def P(Na,Nb,qa,qb):
    return Omega(Na,qa)*Omega(Nb,qb)/Omega(Na+Nb,q)

#Given values
Na = 50
Nb = 50
N = Na+Nb
q  = 100
qa = linspace(0,q,q+1)
#Calculate probabilities and print them in terminal
Prob = zeros(q+1)
for i in qa:
    Prob[i] = P(Na,Nb,i,q-i)
    print('P(N=%s,qa=%s)= '%(N,i),Prob[i])

#Plot probability as function of energy units
plot(qa,Prob)
xlabel('q[units of energy]')
ylabel('P(%s,q)'%N)
title('Probability P(qa) as function of qa for all qa')
show()

max_index, max_value = max(enumerate(Prob), key=operator.itemgetter(1))
print('State with qa = %s most probable with a probability of %s'%(max_index,max_value))
print('State with qb = %s has a probability of %s'%(q-qa[0],Prob[0]))