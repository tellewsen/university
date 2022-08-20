#Exercise 7.6 p. 396

class Quadratic:
#Defines a,b and c in the class.
    def __init__(self,a,b,c):
        self.a = float(a)
        self.b = float(b)
        self.c = float(c)
#Calculates a*x**2 + b*x + c for a chosen value x.
    def value(self,x):
        return self.a*x**2 +self.b*x + self.c

#Prints table of x and y values for n uniformally spaced x values in [L,R]
    def table(self,L,R,n):
        from scitools.std import linspace
        x = linspace(L,R,n)
        y = self.value(x)
        for x,y in zip(x,y):
            print ' %8.4f %8.4f' %(x,y)


#Calculates the roots for a*x**2 + b*x + c (see comment below if curious)
    def roots(self):		
        from numpy import sqrt
        a, b, c = self.a, self.b, self.c
        x1= (-b + sqrt(b**2 - 4*a*c)) / (2*a)
	x2= (-b - sqrt(b**2 - 4*a*c)) / (2*a)
        return x1,x2

quad= Quadratic(1,4,2)
print "Value:", quad.value(0)
quad.table(-2,2,5)
print "Roots:", quad.roots()



#root only works when b**2-4*a*c is positive.
#If you want imaginary numbers there needs to be an if test checking if 
#b**2-4*a*c makes a negative number and then using complex math.
#Unfortunately I have no clue how the cmath module in python works at
#this point...
#There's also the point about numbers being really close to each other and 
#that causing the nominator to be zero even though it shouldn't be.. but 
#you won't run into problems with that unless you really try to do so.



"""
Terminal> python Quadratic.py 
2.0
  -2.0000  -2.0000
  -1.0000  -1.0000
   0.0000   2.0000
   1.0000   7.0000
   2.0000  14.0000
(-0.58578643762690485, -3.4142135623730949)
"""
