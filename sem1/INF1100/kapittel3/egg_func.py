#Exercise 3.15 p.125

from math import pi,log

def egg(M,T0,Ty):
    p= 1.038        #taken from egg.py
    K= 5.4*10**-3   #taken from egg.py
    c= 3.7          #taken from egg.py
    Tw= 100         #taken from egg.py
    return (M**(2./3.)*c*p**(1./3.))/(K*pi**2*(4.*pi/3.)**(2.0/3))*log(0.76*(T0-Tw)/(Ty-Tw))

print 'Hardboiled,Small egg,fridge:', egg(M=47, T0=4, Ty=70), 'Seconds'
print 'Hardboiled,Small egg,room:  ', egg(M=47, T0=25, Ty=70), 'Seconds'
print 'Hardboiled,Large egg,fridge:', egg(M=67, T0=4, Ty=70), 'Seconds'
print 'Hardboiled,Large egg,room:  ', egg(M=67, T0=25, Ty=70), 'Seconds'
print ''
print 'Softboiled,Small egg,fridge:', egg(M=47, T0=4, Ty=63), 'Seconds'
print 'Softboiled,Small egg,room:  ', egg(M=47, T0=25, Ty=63), 'Seconds'
print 'Softboiled,Large egg,fridge:', egg(M=67, T0=4, Ty=63), 'Seconds'
print 'Softboiled,Large egg,room:  ', egg(M=67, T0=25, Ty=63), 'Seconds'

"""
Terminal>  python egg_func.py 
Hardboiled,Small egg,fridge: 313.094549022 Seconds
Hardboiled,Small egg,room:   226.125571496 Seconds
Hardboiled,Large egg,fridge: 396.576342529 Seconds
Hardboiled,Large egg,room:   286.418439338 Seconds

Softboiled,Small egg,fridge: 239.209859774 Seconds
Softboiled,Small egg,room:   152.240882247 Seconds
Softboiled,Large egg,fridge: 302.991449651 Seconds
Softboiled,Large egg,room:   192.83354646 Seconds
"""
