from math import sqrt,pi,exp

x=1.0
m=0.0
s=2.0

f=(1.0/(sqrt(2.0*pi)*s))*exp(-0.5*((x-m)/s)**2.0)

print f



"""
Terminal> python Gaussian1.py
0.176032663382
"""
