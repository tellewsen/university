from random import random

antfeil = 0; N= 10000
x0 = y0 =z0 = 0.0
feildistrib1 = feildistrib2 = 0.0

for apekatt in range(N):
    x = random(); y = random (); z = random()
    distrib1 = (x + y) * z
    distrib2 = x*z + y*z

    if distrib1 != distrib2:
        antfeil +=1
        x0 = x; y0 = y; z0 = z
        feildistrib1 = distrib1
        feildistrib2 = distrib2

print (100. * antfeil/N) #Feilprosenten rett og slett
print (x0, y0, z0, feildistrib1 - feildistrib2) #tilfeldig valgte x,y og z verdier, plus forskjellen mellom (x+y)*z og x*z + y*z for siste utregning
