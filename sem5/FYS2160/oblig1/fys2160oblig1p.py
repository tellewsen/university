# -*- coding: utf-8 -*-
"""
Created on Sat Aug 30 00:27:05 2014

@author: Andreas ELlewsen
"""

import random
import operator
import math
from numpy import zeros,linspace,bincount

#Defines functions
def Omega(N,q):
    return math.factorial(N+q-1)/math.factorial(q)/math.factorial(N-1)

def P(Na,Nb,qa,qb):
    return Omega(Na,qa)*Omega(Nb,qb)/Omega(Na+Nb,q)


#Set constants
N = 120 #Oscillators
q = 120 #Total energy
Na=N/2
Nb=N/2

#Generate initial state
Nlist = zeros(N)
for i in range(q):
    Nlist[random.randint(0,N-1)] +=1

#Nlist[0] = q #Test for all energy at one end

#Develop system in time
n = int(10E5) #number of times to move energy

qa = zeros(n+1,dtype=int)
qb = zeros(n+1,dtype=int)
#time = linspace(0,n,n+1)
i = 0
while i <= n:
    print('%6.4f'%(float(i/n)*100),'% done')
    value1  = random.randint(0,N-1)
    value2  = random.randint(0,N-1)
    if  Nlist[value1] > 0 and value1 != value2:
        Nlist[value1] -= 1
        Nlist[value2] += 1

    qa[i] = sum(Nlist[0:len(Nlist)/2])
#    qb[i] = int(q-qa[i])
    i +=1

#plot(time,qa/Na,time,qb/Nb)
#show()

Prob = bincount(qa)/n
qlist = linspace(0,q,q+1)
Prob = append(Prob,zeros(q+1-len(Prob)))
plot(qlist,Prob)
show()


#Given values
qa = linspace(0,q,q+1)
#Calculate probabilities and print them in terminal
Prob = zeros(q+1)
for i in qa:
    Prob[i] = P(Na,Nb,i,q-i)
    #print('P(N=%s,qa=%s)= '%(N,i),Prob[i])

#Plot probability as function of energy units
plot(qa,Prob)
xlabel('q[units of energy]')
ylabel('P(%s,%s)'%(N,q))
title('Probability P(qa) as function of qa for all qa')
show()

max_index, max_value = max(enumerate(Prob), key=operator.itemgetter(1))
print('State with qa = %s most probable with a probability of %s'%(max_index,max_value))
print('State with qb = %s has a probability of %s'%(q-qa[0],Prob[0]))