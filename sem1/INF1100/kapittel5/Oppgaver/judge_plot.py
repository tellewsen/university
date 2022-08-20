#Exercise 5.18 p. 248

from scitools.std import cos,pi
import numpy as np
x = np.linspace(0,2,1000)
y = cos(18*pi*x)
import matplotlib.pyplot as plt
plt.plot(x,y)
plt.show()


#Here we see that if we only use 20 points, the graph we get looks right but is actually completely wrong. By changing to 1000 points we see a very different graph
