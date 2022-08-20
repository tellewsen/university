from scitools.std import *
import numpy as np
import glob,os
def f(x,t):
    return np.exp(-(x-3*t)**2)*np.sin(3*np.pi*(x-t))



xarray = np.linspace(-6,6,1001)
tarray = np.linspace(-1,1,61)

counter = 0

for t in tarray:
    farray = f(xarray,t)
    plot(xarray,farray,xlabel= 'x',ylabel='f(x,t)',legend='f(x,0)',
         axis=[xarray[0], xarray[-1], -1.1 , 1.1], savefig='tmp%04d.png' %counter)
    counter +=1
movie('tmp*.png')

for filename in glob.glob('tmp*.png'):
    os.remove(filename)
