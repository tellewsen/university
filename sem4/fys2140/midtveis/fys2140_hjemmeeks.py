from numpy import *
from matplotlib.pyplot import *


#We choose the constants from the next problem.
a     = 0.004        	#nm
x0    = 0.127        	#nm
x1    = 0.240        	#nm
omega = 564e3	        #ns-1
M     = 32.5732799552  	#Cl35   mass in eV/c^2
m     = 0.938          	#Proton mass in eV/c^2
my    = (M*m)/(M+m)    	#Reduced mass in eV/c^2

#Arrays
n  = 10E2
x  = linspace(x1 -.05 ,x1+ .05, n)

#Calculate
Psix0 = (1./(2.*pi*a**2.))**(1./2.)*exp(-(x-x1)**2./(2*a**2))
a = 2*a
Psix1 = (1./(2.*pi*a**2.))**(1./2.)*exp(-(x-x1)**2./(2*a**2))
a = 2*a
Psix2 = (1./(2.*pi*a**2.))**(1./2.)*exp(-(x-x1)**2./(2*a**2))

#Plot
hold('on')
plot(x,Psix0,'r-')
plot(x,Psix1,'g-')
plot(x,Psix2,'b-')
title('Sansynlighets tetthet for 3 forskjellige a')
xlabel('Lengde [nm]')
ylabel('Sansynlighetstetthet')
show()
