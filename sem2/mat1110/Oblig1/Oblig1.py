from scitools.std import *

def f1(x,y):
    return 3*y**2 + 2*x*y

def f2(x,y):
    return x**2 - y**2



xarray = linspace(-2,2,101,2)
yarray = linspace(-2,2,101,2)
answer = zeros(len(xarray))
for e in range(len(xarray)):
    for x in xarray:
        for y in yarray:
             answer[e] = 3*y**2 + 2*x*y


plot(xarray,yarray,answer)
