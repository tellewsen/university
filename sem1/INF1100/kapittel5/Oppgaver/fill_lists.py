#Exercise 5.2 p.245

import numpy as np

def h(x):
    return (1.0/np.sqrt(2*np.pi))*np.exp(-0.5*x**2)


a=np.linspace(-4,4,41)

xlist= np.array([i for i in a])
hlist= np.array([h(x) for x in a])

for xlist,hlist in zip(xlist,hlist):
    print '%4.1f  %5f'%(xlist,hlist)


"""
n= 41
dx 1-(n-1)
xlist = np.zeros(41)
hlist = np.zeros(41)
pairs = 
"""
