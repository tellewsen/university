#Exercise 5.17 p. 247
from scitools.std import *
import numpy as np

def f(x,t):
    return np.exp(-(x-3*t)**2)*np.sin(3*np.pi*(x-t))


t=0				#t given in exercise text

xlist = np.linspace(-4,4,250) 	#Makes xlist
flist = np.zeros(len(xlist))	#Makes flist with zeros
for i in range(len(xlist)):	#Fills flist with f(x,t) values
    flist[i] = f(xlist[i],t)

plot(xlist,flist,'-r') 		#Plots the the two arrays
xlabel('x')			
ylabel('f(x,t)')		
legend('f(x,0)')

"""
Terminal> python plot_wavepacket.py 
*plot appears*
"""
