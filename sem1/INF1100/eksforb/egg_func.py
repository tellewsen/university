from math import pi,log
def egg(M,To=20,Ty=70):
    rho = 1.038
    K   = 5.4*10**-3
    c   = 3.7
    Tw  = 100.
    return ((M**(2./3.)*c*rho**(1./3.))/(K*pi**2.*(4*pi/3.)**(2./3.)))*log(0.76*((To-Tw)/(Ty-Tw)))

print 'Hardboiled big egg from fridge  : %f' %(egg(M=67, To=4, Ty=70))
print 'Softboiled big egg from fridge  : %f' %(egg(M=67, To=4, Ty=63))
print 'Hardboiled big egg from room    : %f' %(egg(M=67, To=25, Ty=70))
print 'Softboiled big egg from room    : %f' %(egg(M=67, To=25, Ty=63))

print 'Hardboiled small egg from fridge: %f' %(egg(M=47, To=4, Ty=70))
print 'Softboiled small egg from fridge: %f' %(egg(M=47, To=4, Ty=63))
print 'Hardboiled small egg from room  : %f' %(egg(M=47, To=25, Ty=70))
print 'Softboiled small egg from room  : %f' %(egg(M=47, To=25, Ty=63))


"""
Funker som den skal
"""
