import sys

def f2c(f):
    c= (f-32)*5./9
    return c

try:
    f = float(sys.argv[1])
except IndexError:
    print 'You need to write a temperature in fahrenheit on the command line'
    exit(0)

except ValueError:
    print 'The temperature should be a number'
    exit(0)
print '%.2f degrees fahrenheit is %.2f degrees celsius' %(f, f2c(f))

"""
testet funker
"""
