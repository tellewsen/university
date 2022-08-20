import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

data = np.loadtxt("outputverlet.dat",unpack=True);
xdata = data[0]
ydata = data[1]
zdata = data[2]


#Insert number of planets in data
n_planets = 100
t_final = 5
#Insert number of timesteps to create data
timesteps = int(len(xdata)/n_planets)


#Find indexes of all gravitationally bound objects at last timestep
energydata = np.loadtxt("energies.dat",unpack=True);
kineticdata = energydata[0]
potentialdata = energydata[1]
totaldata = kineticdata + potentialdata/2.
indexes = []
for i in range(0,n_planets):
	if totaldata[i +(timesteps-1)*n_planets] < 0:
		indexes.append(i)

#Empty arrays for positions
x = []
y = []
z = []
masses = np.loadtxt("masses.dat",unpack=True);
M = 0
#Positions of objects arranged in right order
for i in indexes:	
	x.append(xdata[i +(timesteps-1)*n_planets])
	y.append(ydata[i +(timesteps-1)*n_planets])
	z.append(zdata[i +(timesteps-1)*n_planets])
	M +=masses[i]
x  =np.array(x)
y  =np.array(y)
z  =np.array(z)
#Calculate radial distance for all objects
r = np.zeros(np.size(x))
for i in range(np.size(x)):
	r[i] = np.sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i])

#calcualte mass center of cluster

R = 0
for i in range(np.size(r)):
	R += masses[i]*r[i]	#center of mass
R = R/M
for i in range(np.size(r)):
	r[i] = abs(r[i]-R)

#Print average distance from mass center
print 'avg distance = ',np.average(r)
print 'std of r = ', np.std(r)

#Plot number of objects per distance

plt.figure()
plt.hist(r,bins=2*np.arange(0,np.max(r)+1,1))
plt.title('Radial distribution')
plt.ylabel('Number of objects')
plt.xlabel('Radial distance [ly]')
plt.show()

counter = np.zeros((int(np.max(r)))/2+1)
index =0
for j in range(2,int(np.max(r))+1,2):
	for i in range(0,len(r)):
		if j-2< r[i]<j :
			counter[index] +=1
	index +=1

#Calculate number density for each distance
index = 0
for j in range(2,int(np.max(r))+1,2):
	counter[index] = counter[index]/ (4./3*np.pi*((j)**3-(j-2)**3))
	index +=1
j = range(0,int(np.max(r))+1,2)
for i in range(len(j)):
	j[i] = j[i] + 1
j = np.array(j)
print len(indexes)
print j
print len(j)
print counter
print len(counter)


nr = np.linspace(0,np.max(r),200)
def n(nr,n_0,r_0):
	return n_0/(1+(nr/r_0)**4)
n1 = n(nr,0.4,0.5)	
n2 = n(nr,0.8,1.5)
n3 = n(nr,1.2,3.0)
plt.figure()
plt.title('Number density of objects N=%.0f'%(n_planets))
plt.ylabel(r'$n(r)/N^2$')
plt.xlabel(r'$r/N^{-1/3}$')
plt.plot(j/n_planets**(-1./3),counter/n_planets**2,'o')
#plt.plot(nr/n_planets**(-1./3),n1/n_planets**2,'-',label=r'$n_0 = 0.4, r_0 =0.5$')
#plt.plot(nr/n_planets**(-1./3),n2/n_planets**2,'--',label=r'$n_0 = 0.8, r_0 =1.5$')
#plt.plot(nr/n_planets**(-1./3),n3/n_planets**2,'-.',label=r'$n_0 = 1.2,	 r_0 =3.0$')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()
