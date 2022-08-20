#Exercise A.7 p.584
from scitools.std import *

def Fortune(F,p,q,n,I):
    X = zeros(n+1)
    C = zeros(n+1)
    X[0] = F
    C[0] = F*p*q/10.0**4    
    for i in range(1,n+1):
        X[i] = array(X[i-1] + (p/100.0)*X[i-1] - C[i-1])
        C[i] = array(C[i-1] + (I/100.0)*C[i-1])
    return X,C

F = 100. 	# Fortune
p = 5.		# % interest rate
q = 3.		# % to spend
n = 20		# number of years
I = 5.  	# Inflation 

plot(range(n+1), Fortune(F,p,q,n,I)[0])

"""
Terminal>  python fortune_and_inflation1.py 
*plot appears*
"""
