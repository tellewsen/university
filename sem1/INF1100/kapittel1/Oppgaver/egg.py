from math import log,pi

M  = 67                  #Mass big egg
p   = 1.038              #Density in g/cm**3
c   = 3.7                #Specific heat capacity in J/gK
K   = 5.4*10**-3         #Thermal conductivity in 1/10**3 W cm K
Tw  = 100                #Temperature of boiling water
Ty  = 70                 #Temperature when yolk coagulates
T0r = 20                 #Temperature of egg roomtemp
T0f = 4                  #Temperature of egg fridgetemp 


tr = (M**(2./3.)*c*p**(1./3.))/(K*pi**2*(4.*pi/3.)**(2.0/3))*log(0.76*(T0r-Tw)/(Ty-Tw))
tf = (M**(2.0/3)*c*p**(1.0/3))/(K*pi**2*(4.*pi/3.)**(2.0/3))*log(0.76*(T0f-Tw)/(Ty-Tw))
print "Time yolk reach %g Celsius if egg roomtemp  : %g seconds" %(Ty,tr)
print "Time yolk reach %g Celsius if egg fridgetemp: %g seconds" %(Ty,tf)

"""
Terminal> $ python egg.py
Time yolk reach 70 Celsius if egg roomtemp  : 315.218 seconds
Time yolk reach 70 Celsius if egg fridgetemp: 396.576 seconds
"""
