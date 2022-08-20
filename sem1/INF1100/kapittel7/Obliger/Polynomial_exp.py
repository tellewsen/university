#Exercise 7.29 p. 404
import Polynomial, sys
from math import exp,factorial

"Remember to enter x as first argument and any number of Ns after that"
"Example: Terminal> python Polynomial_exp.py x N N N N...."


x = float(sys.argv[1])
Na = sys.argv[2:]
N= []
for i in Na:
    N.append( int(i) )


print 'Exact : %.16f' %(exp(x))
print ' N  Approx'

for i in N:
    px = []; factor = 1
    for k in range(1,i+1):
        px.append(1./factor)
        factor *= k
    poly = Polynomial.Polynomial(px)
    print '%2g %20.16f'%(i,poly(x))


""""
Terminal> python Polynomial_exp.py 0.5 2 5 10 15 25
Exact : 1.6487212707001282
 N  Approx
 2   1.6250000000000000
 5   1.6486979166666667
10   1.6487212706873655
15   1.6487212707001278
25   1.6487212707001278

Terminal>  python Polynomial_exp.py 3 2 5 10 15 25
Exact : 20.0855369231876679
 N  Approx
 2   8.5000000000000000
 5  18.3999999999999986
10  20.0796651785714246
15  20.0855344309708101
25  20.0855369231876537

Terminal> python Polynomial_exp.py 10 2 5 10 15 25
Exact : 22026.4657948067178950
 N  Approx
 2  61.0000000000000000
 5 1477.6666666666665151
10 12842.3051146384477761
15 20952.8869686065445421
25 22026.0763608910674520
"""
