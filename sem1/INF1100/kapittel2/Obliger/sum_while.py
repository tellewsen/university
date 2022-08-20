s = 0; k = 1; M = 100

while k <= M:    #Mistake number 1, did not include M, fixed by changing from < to <=
    s += 1./k    #Mistake number 2, integer division, fixed by adding a .
    k += 1       #Mistake number 3, did not increase k so it never reached M. Made it go forever.
print s


"""
Kjoreeksempel>  python sum_while.py 
5.18737751764
"""
