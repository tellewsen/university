#Exercise 8.6 p. 467

import random,sys

"""
Computes the probability of getting at least one 6 when throwing n dice N times
Usage: "python one6_ndice.py n N" (where n is number of dice and N is number of times to do the experiment)
"""
K=0


n = int(sys.argv[1])	#Number of dice
N = int(sys.argv[2])	#Times done

for i in range(N):	#Does experiment N times
    r = 0
    M = 0
    for i in range(n):	#Throws dice n times
        r = (random.randint(1,6))
        if r == 6:
            M+=1
    if M > 0:
        K += 1

prob = float(K)/N

print 'Numerical prob: %f' %(prob)
print 'Exact prob     :', 11./36, '(ignore this line if not testing for 2 dice)'

"""
Terminal> python one6_ndice.py 2 1000
Numerical prob: 0.280000
Exact prob     : 0.305555555556 (ignore this line if not testing for 2 dice)

Terminal> python one6_ndice.py 2 1000
Numerical prob: 0.304000
Exact prob     : 0.305555555556 (ignore this line if not testing for 2 dice)

Terminal> python one6_ndice.py 2 100000
Numerical prob: 0.306850
Exact prob     : 0.305555555556 (ignore this line if not testing for 2 dice)

Terminal> python one6_ndice.py 4 1000
Numerical prob: 0.488000
Exact prob     : 0.305555555556 (ignore this line if not testing for 2 dice)
"""
