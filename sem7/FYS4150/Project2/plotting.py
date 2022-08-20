from matplotlib.pyplot import *
from numpy import *
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

#Part 1

data = loadtxt("output.txt",unpack=True);

figure()
plot(data[0],data[1])
title("Energies for two electrons in potential well")
xlabel(r'$\rho$')
ylabel(r'E($\rho$)')
show()

#Part 2 
"""
data1 = loadtxt("output1.txt",unpack=True);
data2 = loadtxt("output2.txt",unpack=True);
data3 = loadtxt("output3.txt",unpack=True);
data4 = loadtxt("output4.txt",unpack=True);

figure(0)
plot(data1[0],(data1[1]))
hold('on')
plot(data1[0],(data2[1]))
plot(data1[0],(data3[1]))
plot(data1[0],(data4[1]))

legend([r"$\omega = 0.01$",r"$\omega = 0.5$",r"$\omega = 1$",r"$\omega = 5$"],loc='best')
title(r'Wave function $\psi(\rho)$')
xlabel(r'$\rho$')
ylabel(r'$\psi(\rho)$')
hold('off')
"""
#Part 3
"""
data1 = loadtxt("output1.txt",unpack=True);
data2 = loadtxt("output2.txt",unpack=True);
data3 = loadtxt("output3.txt",unpack=True);
data4 = loadtxt("output4.txt",unpack=True);

figure(1)
plot(data1[0],(data1[1])**2)
hold('on')
plot(data1[0],(data2[1])**2)
plot(data1[0],(data3[1])**2)
plot(data1[0],(data4[1])**2)

legend([r"$\omega = 0.01$",r"$\omega = 0.5$",r"$\omega = 1$",r"$\omega = 5$"],loc='best')
title(r'Probability distribution $|\psi(\rho)|^2$')
xlabel(r'$\rho$')
ylabel(r'$|\psi(\rho)|^2$')
show()
"""
