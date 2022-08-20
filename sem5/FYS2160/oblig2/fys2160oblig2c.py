# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 15:21:59 2014

@author: Andreas Ellewsen
"""

import random

N = 50   #Oscillators
n = 10000   #Iterations
Nlist = np.zeros((n,N),dtype=int)
Macrostate = np.zeros(n)

#Generate n microstates with N Oscillators
for i in range(n):
    for j in range(N):
        Nlist[(i,j)] = random.choice([-1,1])
    Macrostate[i] = sum(Nlist[(i)])

s = linspace(-30,30)
Omega = exp(-2*s**2/N)
hist(Macrostate,26)
xlabel('Total Energy/[-µB] ')
ylabel('Number of trial')
plot(s,Omega)
show()