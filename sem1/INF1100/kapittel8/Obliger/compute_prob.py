#Exercise 8.2 p. 466

import random,sys

for i in 1,2,3,6:	#N = 1,2,3,,6
    N = 10**i		#N = 10,100,1000,100000
    r = []
    M = float(0)

    for i in range(N):	#fills r with N values in [0,1)
        r.append(random.random())

    for i in range(len(r)):	#Runs through r and increases M by 1 if 0.5 < r[i] < 0.6
       if 0.5 < r[i] < 0.6:
          M +=1
    print ''
    print 'Drew %g times, was in (0.5,0.6) %g times' %(N,M)	#Prints N and M
    print 'Probability(100*M/N) = %f%%' %(100*M/N)		#Prints probability(100*M/N)

"""
Terminal> python compute_prob.py

Drew 10 times, was in (0.5,0.6) 3 times
Probability(100*M/N) = 30.000000%

Drew 100 times, was in (0.5,0.6) 16 times
Probability(100*M/N) = 16.000000%

Drew 1000 times, was in (0.5,0.6) 90 times
Probability(100*M/N) = 9.000000%

Drew 1e+06 times, was in (0.5,0.6) 99888 times
Probability(100*M/N) = 9.988800%
"""
