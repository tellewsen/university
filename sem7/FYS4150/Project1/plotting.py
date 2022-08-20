from matplotlib.pyplot import *
from numpy import *

data = loadtxt("output.txt",unpack=True);

figure()
plot(data[0],data[1])
hold("on")
plot(data[0],data[2])
plot(data[0],data[3])
plot(data[0],data[4])
legend(["u(x)","v_special(x)","v_gauss(x)","v_LU(x)"])
title("Solutions of the Poisson Equation")
xlabel('x')
ylabel('y(x)')
show()
