#Exercise 6.19 p. 335

import sys
from scitools.std import *

f = StringFunction(str(sys.argv[1]))
f.vectorize(globals())
a = float(eval(sys.argv[2]))
b = float(eval(sys.argv[3]))
n = int(sys.argv[4])
output = open(sys.argv[5],'w')

x = linspace(a,b,n)
y = f(x)

for x,y in zip(x,y):
    output.write('%f %f \n' %(x,y))

output.close()

"""
Terminal> python write_cml_function.py 'cos(x)' 0 2*pi 5 'output.dat'
"""
