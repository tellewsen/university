#Exercise 9.11 p. 551
from numpy import *
from scitools.std import *
from Diff import *

class Backward2(Diff):
    def __call__(self,x):
        f,h = self.f, self.h
        return (f(x-2*h) - 4*f(x-h) + 3*f(x))/(2*h)

#test below



if __name__ == '__main__':

    def f(t):
        return exp(-t)

    def exact(t):
        return -exp(-t)


    for k in range(15):
        h = 2**-k
        b1 = Backward1(f,h)
        b2 = Backward2(f,h)
        print 'k=%2g          h= %g' %(k,h)
        print 'Backward1 error: %.2e'%(abs(exact(0) - b1(0)))
        print 'Backward2 error: %.2e'%(abs(exact(0) - b2(0)))
        print ''
"""
Terminal> python Backward2.py 
k= 0          h= 1
Backward1 error: 7.18e-01
Backward2 error: 7.58e-01

k= 1          h= 0.5
Backward1 error: 2.97e-01
Backward2 error: 1.23e-01

k= 2          h= 0.25
Backward1 error: 1.36e-01
Backward2 error: 2.52e-02

k= 3          h= 0.125
Backward1 error: 6.52e-02
Backward2 error: 5.73e-03

k= 4          h= 0.0625
Backward1 error: 3.19e-02
Backward2 error: 1.36e-03

k= 5          h= 0.03125
Backward1 error: 1.58e-02
Backward2 error: 3.33e-04

k= 6          h= 0.015625
Backward1 error: 7.85e-03
Backward2 error: 8.23e-05

k= 7          h= 0.0078125
Backward1 error: 3.92e-03
Backward2 error: 2.05e-05

k= 8          h= 0.00390625
Backward1 error: 1.96e-03
Backward2 error: 5.10e-06

k= 9          h= 0.00195312
Backward1 error: 9.77e-04
Backward2 error: 1.27e-06

k=10          h= 0.000976562
Backward1 error: 4.88e-04
Backward2 error: 3.18e-07

k=11          h= 0.000488281
Backward1 error: 2.44e-04
Backward2 error: 7.95e-08

k=12          h= 0.000244141
Backward1 error: 1.22e-04
Backward2 error: 1.99e-08

k=13          h= 0.00012207
Backward1 error: 6.10e-05
Backward2 error: 4.97e-09

k=14          h= 6.10352e-05
Backward1 error: 3.05e-05
Backward2 error: 1.24e-09

"""
