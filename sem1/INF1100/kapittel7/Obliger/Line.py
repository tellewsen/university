#Exercise 7.4 p. 395

class Line:
    def __init__(self,p1,p2):
        self.a = (p2[1] -p1[1]) / float(p2[0] - p1[0])
        self.b = p1[1] - self.a*p1[0]

    def value(self,x):
        return self.a*x + self.b

#Test from book:
line1 = Line((0,-1),(2,4))
print line1.value(0.5), line1.value(0), line1.value(1)


#Test, easy to see that it works:
line2 = Line((0,0), (1,2))
print line2.value(1), line2.value(2), line2.value(3)



"""
Terminal> python Line.py 
0.25 -1.0 1.5
2.0 4.0 6.0
"""
