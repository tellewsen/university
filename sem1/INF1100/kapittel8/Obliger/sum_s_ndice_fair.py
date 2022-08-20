#Exercise 8.9 p 467

import random,sys


#Try/except block for testing that you input N and r. 
#Describes usage if you fail to input values.
try:
    N = int(sys.argv[1])		#Times to do experiment		
    q = float(sys.argv[2])		#Price of playing
    s = float(sys.argv[3])		#Sum die should be less than
except IndexError:
    print "Usage: python sum9_s_ndice_fair.py N q s \
          (N=times to do experiment, q=cost of playing, \
          s=what the number of eyes should be less than)"
    exit(0)

win  = 0				#Win counter

for i in range(N):
    dicesum = 0				#Counter sum of dice
    for i in range(1,5):		#Throws dice 4 times
        dicesum += random.randint(1,6)	#Throws dice
    if dicesum < s:			#If sum of the 4 throws is less than 9
        win +=1				#increases money by 9
    
try :
    p = float(win)/N			#Calculates win percentage
    r = float(q)/p			#Calculates what r is fair
except ZeroDivisionError:
    print 'Won 0 times, unable to compute r, please increase N and/or s'


try :
    print 'By using the values you inserted you would get a fair game if r was %s' %(r)
except NameError:
    exit(0)

"""
Terminal> python sum_s_ndice_fair.py 1000000 1 9
By using the values you inserted you would get a fair game if r was 18.5253797703

Terminal> python sum_s_ndice_fair.py 100000 1 9
By using the values you inserted you would get a fair game if r was 18.8359389716

Terminal> python sum_s_ndice_fair.py 1000 1 9
By using the values you inserted you would get a fair game if r was 21.7391304348

Terminal> python sum_s_ndice_fair.py 100 1 9
By using the values you inserted you would get a fair game if r was 33.3333333333

Terminal> python sum_s_ndice_fair.py 10 1 9
Won 0 times, unable to compute r, please increase N and/or s

Terminal> python sum_s_ndice_fair.py 10 1  
Usage: python sum9_s_ndice_fair.py N q s           (N=times to do experiment, q=cost of playing,           s=what the number of eyes should be less than)
"""
