#Exercise 9.2 p. 549

from Parabola import *
class Cubic(Parabola):
    def __init__(self,c0,c1,c2,c3):
        Parabola.__init__(self,c0,c1,c2)
        self.c3 = c3

    def __call__(self,x):
        print 'Call to Cubic'
        return Parabola.__call__(self,x) + self.c3*x**3


class Poly4(Cubic):
    def __init__(self,c0,c1,c2,c3,c4):
        Cubic.__init__(self,c0,c1,c2,c3)
        self.c4 = c4

    def __call__(self,x):
        print 'Call to Poly4'
        return Cubic.__call__(self,x) + self.c4*x**4


#Test below
if __name__ == "__main__":

    print 'Use of Cubic starts here'
    a = Cubic(1,2,3,4)
    print a(1)
    print 'Calling table function from Cubic here'
    print a.table(1,2,5)

    print 'Use of Poly4 starts here'
    b= Poly4(1,2,3,4,5)
    print b(1)
    print 'Calling table function from Poly4 here'
    print b.table(1,2,5)

"""
Terminal> python Cubic_Poly4.py 
Use of Cubic starts here
Call to Cubic
10
Calling table function from Cubic here
Call to Cubic
Call to Cubic
Call to Cubic
Call to Cubic
Call to Cubic
           1           10
        1.25           16
         1.5        24.25
        1.75       35.125
           2           49

Use of Poly4 starts here
Call to Poly4
Call to Cubic
15
Calling table function from Poly4 here
Call to Poly4
Call to Cubic
Call to Poly4
Call to Cubic
Call to Poly4
Call to Cubic
Call to Poly4
Call to Cubic
Call to Poly4
Call to Cubic
           1           15
        1.25       28.207
         1.5      49.5625
        1.75      82.0195
           2          129

"""
