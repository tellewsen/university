#Exercise 3.24 p.128

def H(x):		#Defines function H(x)
    if x<0:		
        return 0	#If x is less than 0 returns 0
    else:
        return 1	#If x is 0 or higher returns 1

print 'H(-0.5) = %g' %(H(-0.5))	#Prints the return from the function when x= -0.5 (Should get 0)
print 'H(0)    = %g' %(H(0))	#Prints the return from the function when x=  0 (Should get 1)
print 'H(10)   = %g' %(H(10))	#Prints the return from the function when x= 10 (Should get 1)


"""
Terminal> python Heaviside.py 
H(-0.5) = 0
H(0)    = 1
H(10)   = 1
"""
