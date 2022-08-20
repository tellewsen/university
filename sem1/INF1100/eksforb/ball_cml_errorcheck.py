import sys

v0 = float(sys.argv[1]) 
g  = float(sys.argv[2]) 
t  = float(sys.argv[3]) 

if 0 < t < 2*v0/g:
    y  = v0*t - 0.5*g*t**2
    print y
else:
    print 't was not in the interval (0,2*v0*g)'

"""
testet funker
"""
