import sys

def f2c(f):
    c= (f-32)*5./9
    return c

f = float(sys.argv[1])
print '%.2f degrees fahrenheit is %.2f degrees celsius' %(f, f2c(f))

"""
testet funker
"""
