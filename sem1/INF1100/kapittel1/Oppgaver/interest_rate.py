A = 1000 #euros
p = 5.0 #5% interest
n = 3.0 #years
interest = A*((1.0+p/100.0)**n)
print "%g euros with %g %% interest will grow to %g euros after %g years" %(A,p,interest,n)

"""
Samplerun> python interest_rate.py 
1000 euros with 5 % interest will grow to 1157.63 euros after 3 years
1157.625
"""
