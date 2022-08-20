#Exercise 5.2 p.244

from math import sqrt,pi,exp

def h(x):						#The function given in the text
    return (1.0/sqrt(2*pi))*exp(-0.5*x**2)

xlist = [i*.2 for i in range(-20,21)]			#Makes a list of 41 equally spaced elements in the interval [-4,4]
hlist = [h(x) for x in xlist]				#List of h(x) values for all elements in xlist

print '   x      h(x)'
for x,h in zip(xlist, hlist):				#Prints so you get something to look at.
    print '%5.1f  %6.7f'   %(x,h)


"""
Terminal>  python fill_arrays_loop.py 
   x      h(x)
 -4.0  0.0001338
 -3.8  0.0002919
 -3.6  0.0006119
 -3.4  0.0012322
 -3.2  0.0023841
 -3.0  0.0044318
 -2.8  0.0079155
 -2.6  0.0135830
 -2.4  0.0223945
 -2.2  0.0354746
 -2.0  0.0539910
 -1.8  0.0789502
 -1.6  0.1109208
 -1.4  0.1497275
 -1.2  0.1941861
 -1.0  0.2419707
 -0.8  0.2896916
 -0.6  0.3332246
 -0.4  0.3682701
 -0.2  0.3910427
  0.0  0.3989423
  0.2  0.3910427
  0.4  0.3682701
  0.6  0.3332246
  0.8  0.2896916
  1.0  0.2419707
  1.2  0.1941861
  1.4  0.1497275
  1.6  0.1109208
  1.8  0.0789502
  2.0  0.0539910
  2.2  0.0354746
  2.4  0.0223945
  2.6  0.0135830
  2.8  0.0079155
  3.0  0.0044318
  3.2  0.0023841
  3.4  0.0012322
  3.6  0.0006119
  3.8  0.0002919
  4.0  0.0001338
"""
