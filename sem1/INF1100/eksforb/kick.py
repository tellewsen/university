from math import pi


rho = 1.2 		#kg/m3	Density of air
V   = [10/3.6,120/3.6] 	#m/s	Velocity of object
a   = 0.11		#m	Radius of football
A   = pi*a**2 		#m2	Cross sectional area
Cd  = 0.2 		#	Drag coeffisient
m   = 0.43		#	Mass of football
g   = 9.81		#m/s2	Velocity of gravity
Fg  = m*g		#N	Gravity force on object

Fd = [0,0]
for i in 0,1:
    Fd[i]  = 0.5*Cd*rho*A*V[i]**2	#	Dragforce due to air on object

print "Kick       DragForce   Ratio Force vs Gravity"
print "Softkick:  %.1f         %.1f" %(Fd[0], Fd[0]/Fg)
print "Hardkick:  %.1f         %.1f" %(Fd[1], Fd[1]/Fg)
