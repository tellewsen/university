from math import pi

Cd=0.2          #Drag coefficient
q=1.2           #Density of air in       kg/m**3
a=0.11          #Radius of football in   m
A=pi*a**2.0     #Cross-sectional area    m**2
Vsoft=120/3.6   #Velocity in             m/s
Vhard=10/3.6    #Velocity in             m/s
g=9.81          #Gravity in              m/s**2
m=0.43          #Mass of football in     kg

Fg=m*g                        #Force of gravity
Fdhard=0.5*Cd*q*A*Vsoft**2.0  #Force of Hardkick
Fdsoft=0.5*Cd*q*A*Vhard**2.0  #Force of Softkick
Ratiohard= Fdhard/Fg          #Ratio of Hardkick/Gravity
Ratiosoft= Fdsoft/Fg          #Ratio of Softkick/Gravity

print """Softkick:
Force of Softkick(10km/h) = %.1f Newton
Force of Gravity  = %.1f Newton
Ratio Softkick/Gravity= %.1f
""" %(Fdsoft,Fg,Ratiosoft)

print """Hardkick:
Force of Hardkick(120km/h) = %.1f Newton
Force of Gravity  = %.1f Newton
Ratio Hardkick/Gravity= %.1f
""" %(Fdhard,Fg,Ratiohard)


"""
Samplerun> python kick.py 
Softkick:
Force of Softkick(10km/h) = 0.0 Newton
Force of Gravity  = 4.2 Newton
Ratio Softkick/Gravity= 0.0

Hardkick:
Force of Hardkick(120km/h) = 5.1 Newton
Force of Gravity  = 4.2 Newton
Ratio Hardkick/Gravity= 1.2
"""
