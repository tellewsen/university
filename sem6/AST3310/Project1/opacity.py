from numpy import *

#Function that reads the kappa table
lines = genfromtxt('opacity.txt')


T    = 1E5
rho  = 1E3
rho  = rho*1E3	 #Convert from kg/m**3 to g/cm**3
R    = rho/(T/1E6) #Convert R to the form used in opacity.txt

Rvalue = log10(R) 
print Rvalue
Tvalue = log10(T)
print Tvalue

def round_to(n, precision):
    correction = 0.5 if n >= 0 else -0.5
    return int( n/precision+correction ) * precision



#Round the both values to the closest value found in the table
Rvalue = (round(2*Rvalue)/2.)

if Tvalue <= 6:
	Tvalue = round_to(Tvalue, 0.05)

if Tvalue > 6:
	Tvalue = round_to(Tvalue, 0.1)
 
 
print Rvalue
print Tvalue
#Picks the right R value from the table 
for i in range(len(lines[0])):
       if Rvalue == lines[0][i]:
        Ri = i        
        print 'Ri:',Ri

#Pick the right T value from the table
for i in range(len(lines)):
    if Tvalue ==lines[i][0]:
        Ti = i
        print 'Ti:',Ti
    print i


#Pick the right kappa value with R and T
kappa = lines[Ti][Ri]
print Tvalue
print lines[70][0]
print 'kappa=',kappa