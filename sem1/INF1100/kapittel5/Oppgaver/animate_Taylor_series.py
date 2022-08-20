#Exercise 5.28 p. 251
from sys import *
from scitools.std import *
import numpy as np
def animate_series(fk, M, N, xmin, xmax, ymin, ymax, n, exact):
    return fk

xarray=np.linspace(sys.argv[4],sys.argv[5],sys.argv[8])

yarray=np.zeros(len(xarray))
yarray=animate_series(xarray)

counter = 0
for i in range(M,N):
    plot(xarray,yarray,'b-',
        xlabel = 'x',
        ylabel = 'y',
        legend = functionhere,
        axis   = [sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7]],
        title  = fk,
        savefig= 'tmp%04d.png' %counter)
    counter += 1
#movie('tmp*.png')
