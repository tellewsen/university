#Exercise A.18 p.
from scitools.std import *

def cosine_Taylor(x,n):
    sp = 0     
    ap = 1
    for j in range(1,n+2):
        s = sp + ap
        a = ap * (-1.)*x**2 / ((2*j-1)*(2*j))
        ap = a
        sp = s
    return s,abs(a)

print 'test for x*pi/2 with x=0->20 and n=25:'
for i in range(0,21):
    answer = cosine_Taylor(i*pi/2,25)[0]
    print '%2.f*pi/2   %20.17f' %(i,answer)

print 'test for pi/2 with n=0 -> 20:'
for i in range(0,21):
    answer = cosine_Taylor(pi/2,i)[0]
    print '   n=%2.f   %20.17f' %(i,answer)
