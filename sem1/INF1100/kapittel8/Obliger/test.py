class velocity:
    def __init__(self,g,v0):
        self.g  = g
        self.v0 = v0
    def __call__(self,t):
        #g,v0 = self.g,self.v0
        return 1./2*self.g*self.v0**2*t
    


class apekatt(velocity):
    def __init__(self,g,v0,r):
       velocity.__init__(self,g,v0)
       self.r = r
    def __call__(self,t):
        return velocity.__call__(self,t) + self.r

a = velocity(9.81,1)
print a(2)

b = velocity(2,2)
print b(4)

c = apekatt(3,4,5)
print c(6)
