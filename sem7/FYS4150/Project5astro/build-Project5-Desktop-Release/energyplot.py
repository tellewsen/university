import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing




energydata = np.loadtxt("energies.dat",unpack=True);
kineticdata = energydata[0]
potentialdata = energydata[1]
totaldata = kineticdata + potentialdata/2
#Insert number of planets in data
n_planets = 100
t_final = 10
#Insert number of timesteps to create data
timesteps = int(len(kineticdata)/n_planets)



kinetic = np.zeros(timesteps)
potential = np.zeros(timesteps)
total = np.zeros(timesteps)
for t in range(0,timesteps):
	for i in range(0,n_planets):	
		kinetic[t] += kineticdata[i +t*n_planets]
		potential[t] += potentialdata[i +t*n_planets]/2
	total[t] += kinetic[t]+potential[t]
plt.figure()
t =  np.linspace(0,t_final,timesteps)
plt.plot(t,kinetic/abs(total[0]),'-')
plt.plot(t,potential/abs(total[0]),'--')
plt.plot(t,total/abs(total[0]),'-.')
plt.xlim([0,t_final])
plt.title('Energy of the system all objects included')
plt.legend(['kinetic','potential','total'],loc='best')
plt.xlabel(r'time[$\tau_{crunch}$]')
plt.ylabel(r'Energy/$|E_0|$')

#Find indexes of all gravitationally bound objects at last timestep
indexes = []
for i in range(0,n_planets):
	if totaldata[i +(timesteps-1)*n_planets] < 0:
		indexes.append(i)
kinetic2 = np.zeros(timesteps)
potential2 = np.zeros(timesteps)
total2 = np.zeros(timesteps)
for t in range(0,timesteps):
	for i in indexes:
		kinetic2[t] += kineticdata[i +t*n_planets]
		potential2[t] += potentialdata[i +t*n_planets]/2
	total2[t] = kinetic2[t]+potential2[t]
plt.figure()
t =  np.linspace(0,t_final,timesteps)
plt.plot(t,kinetic2/abs(total2[0]),'-')
plt.plot(t,potential2/abs(total2[0]),'--')
plt.plot(t,total2/abs(total2[0]),'-.')
plt.xlim([0,t_final])
plt.title('Energy of the system(only grav. bound)')
plt.legend(['kinetic','potential','total'],loc='best')
plt.xlabel(r'time[$\tau_{crunch}$]')
plt.ylabel(r'Energy/$|E_0|$')
plt.show()

print total2[-1]-total2[0]
print abs((total2[-1]-total2[0]))/abs(total2[0])

print "2* avg kinetic = ",2*kinetic2[-1]
print "-avg potential = ",-potential2[-1]

#Test of virial theorem
plt.figure()
plt.plot(t,(2*kinetic2+potential2)/abs(total2[0]))
plt.title('Virial theorem')
plt.xlabel(r'time[$\tau_{crunch}$]')
plt.ylabel(r'$2<K>/|E_0| + <V>/|E_0|$')
plt.show()

#Calculate time average instead

#average kinetic energy
avgkin = 0
avgpot = 0
for t in xrange(int(timesteps*3./4),timesteps):
	avgkin += kinetic2[t]
	avgpot += potential2[t]

avgkin = avgkin/(timesteps-int(timesteps*3./4))
avgpot = avgpot/(timesteps-int(timesteps*3./4))
print 2*avgkin
print -avgpot
