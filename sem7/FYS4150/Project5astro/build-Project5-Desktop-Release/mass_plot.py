import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

"""
This function plots the mass distribution of the system.
"""

data = np.loadtxt("masses.dat",unpack=True);
plt.figure(0)
plt.hist(data,bins=20,normed=1)
plt.title(r'Mass distribution')
plt.ylabel(r'P(M)',size= 14)
plt.xlabel(r'Mass [kg]',size=14)
plt.legend(loc='best')
plt.xlim([np.min(data),np.max(data)])
plt.show()

print "M_tot = ",sum(data)
