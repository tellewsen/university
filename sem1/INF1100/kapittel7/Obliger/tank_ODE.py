#Exercise E.3 p. 707

R	#Radius tank
h0 	#height
r	#Radius valve
h(t)	#decreseases with time

in time dt reduce by dh
hight corresponds to volume pi*R**2*dh
water corresponds to pi*r**2*v*dt	#v= outflow velocity


def v(t):
    return sqrt(2*g* h(t) + hd(t)**2)
