#Exercise A.5 p.586

x = 100				# initial amount
p = 5				# interest rate
n = 4				# number of years
count = 0

xp= x				# previous x

print 'years amount'
while count <= n:
    x  = xp + (p/100.0)*xp
    count +=1
    xp = x
    print '%5.0f %6f'%(count,x)


"""
Terminal>  python growth_years_efficient.py 
years amount
    1 105.000000
    2 110.250000
    3 115.762500
    4 121.550625
    5 127.628156
"""
