#Exercise 8.8 p.467

import random,sys


#Try/except block for testing that you input N and r. 
#Describes usage if you fail to input values.
try:
    N = int(sys.argv[1])			#Times to do experiment		
    r = float(sys.argv[2])			#Prize for getting less than 9 in 4 throws

except IndexError:
    print "Usage: python sum9_4dice.py N r (N=times to do experiment, r= prize for winning)"
    exit(0)

except ValueError:
    print "N must be an integer"
    exit(0)

money = 0				#Startmoney in euros

for i in range(N):
    dicesum = 0				#Counter sum of dice
    for i in range(1,5):		#Throws dice 4 times
        dicesum += random.randint(1,6)	#Throws dice
    if dicesum < 9:			#If sum of the 4 throws is less than 9
        money +=r-1			#increases money by 9
    else:				#If sum is 9 or more decreases by 1
        money -=1


if money < 0 :
    print "Threw %g times, ended up with %g euro." %(N,money)
    print "This is a good way to lose all your money!"

elif money > 0 :
    print "Threw %g times, ended up with %g euro." %(N,money)
    print "You're gonna get RICH!"

else:
    print "Threw %g times, ended up with %g euro." %(N,money)
    print "Looks like you end up where you started"

"""
Terminal> sum9_4dice.py 10 10
Threw 10 times, ended up with 10 euro.
You're gonna get RICH!

Terminal> python sum9_4dice.py 100 10
Threw 100 times, ended up with -10 euro.
This is a good way to lose all your money!

Terminal> python sum9_4dice.py 1000 10
Threw 1000 times, ended up with -390 euro.
This is a good way to lose all your money!

Terminal> python sum9_4dice.py 10000 10
Threw 10000 times, ended up with -4130 euro.
This is a good way to lose all your money!

Terminal> python sum9_4dice.py 100     
Usage: python sum9_4dice.py N r (N=times to do experiment, r= prize for winning)

Terminal> python sum9_4dice.py 10.2 10
N must be an integer
"""
