import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

"""
This is the plotting file used for output using method 1.
Note the "output.dat" on the line below indicating the datafile used.
You MUST set this yourself.
This file is also set up so that you get a plot with the percent of accepted moves,
and because of this the user MUST update the variables mcs and spins with the values used for the datafile in use.
"""

mcs = 1000000.
spins = 20.
totalmoves = spins*spins*mcs


data = np.loadtxt("output.dat",unpack=True);
plt.figure(0)
plt.plot(data[0],data[1]/(totalmoves)*100)
plt.title(r'Percentage of moves accepted for each temperature')
plt.ylabel(r'Percent accepted moves',size= 14)
plt.xlabel(r'Temperature [kT/J]',size=14)
plt.legend(loc='best')
plt.xlim([np.min(data[0]),np.max(data[0])])
plt.show()
