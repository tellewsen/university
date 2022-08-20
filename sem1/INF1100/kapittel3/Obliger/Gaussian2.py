#Exercise 3.14 p. 124

from math import sqrt,pi,exp  #Imports functions from math needed for function gauss

def gauss(x,m=0,s=1):                                 #Defines function gauss
    f=(1.0/(sqrt(2.0*pi)*s))*exp(-0.5*((x-m)/s)**2.0) #The equation gauss uses (f)
    return f                  #Returns answer for f

x=-5                          #Startvalue for x
while x<=5:                   #Checks that x has not exceeded 5 and stops it when it does
    print '%2.0f = %g' %(x,gauss(x))            #Prints answers to gauss(x)
    x += 1                    #Increases x by 1 before loop starts again.



"""
Terminal> python Gaussian2.py 
-5 = 1.48672e-06
-4 = 0.00013383
-3 = 0.00443185
-2 = 0.053991
-1 = 0.241971
 0 = 0.398942
 1 = 0.241971
 2 = 0.053991
 3 = 0.00443185
 4 = 0.00013383
 5 = 1.48672e-06

Vet ikke helt om oppgaven sier at jeg skal skrive ut for -5, -4 ,3 .., 5 eller -5, -4.9, -4.8, ..., 4.9, 5.0
Har satt opp for at den skal bruke alternativ 1, men om det er alternativ 2 som er riktig endres siste linje i koden til "x +=.1"
og nest siste linje til "print '%2.2f = %g' %(x, gauss(x))"
"""
