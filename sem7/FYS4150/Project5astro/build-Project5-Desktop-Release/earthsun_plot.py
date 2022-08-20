#from scitools.easyviz import movie
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

"""
This file makes a video of the position of all objects for each time step
"""
data = np.loadtxt("outputRK4.dat",unpack=True);
data2 = np.loadtxt("outputverlet.dat",unpack=True);
xdata = data[0]
ydata = data[1]
zdata = data[2]
xdata2 = data2[0]
ydata2 = data2[1]
zdata2 = data2[2]
#Insert number of planets in data
n_planets = 2

#Insert number of timesteps to create data
timesteps = int(len(xdata)/n_planets)

AU = 1.496E11
ly = 9.4607e15

x = np.zeros((n_planets,timesteps))
y = np.zeros((n_planets,timesteps))
z = np.zeros((n_planets,timesteps))
x2 = np.zeros((n_planets,timesteps))
y2 = np.zeros((n_planets,timesteps))
z2 = np.zeros((n_planets,timesteps))

for t in range(0,timesteps):
	for i in range(0,n_planets):	
		x[i,t] = xdata[i +t*n_planets]/AU
		y[i,t] = ydata[i +t*n_planets]/AU
		z[i,t] = zdata[i +t*n_planets]/AU
		x2[i,t] = xdata2[i +t*n_planets]/AU
		y2[i,t] = ydata2[i +t*n_planets]/AU
		z2[i,t] = zdata2[i +t*n_planets]/AU

xmax = 1.2
ymax = 1.2
zmax = 1.2

plt.figure()
plt.plot(x[1,:],y[1,:])
plt.plot(x2[1,:],y2[1,:],'--')
#plt.title('Velocity-Verlet')
plt.title('Comparison of RK4 and Velocity Verlet')
plt.legend(['RK','V-V'])
plt.xlim([-xmax,xmax])
plt.ylim([-ymax,ymax])
plt.xlabel('X[AU]')
plt.ylabel('Y[AU]')
plt.show()


"""
#Plot the trajectories of the objects
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
for i in range(0,n_planets):
	ax.plot(x[i,:],y[i,:],z[i,:])
plt.title('Velocity-Verlet')
#plt.title('Runge Kutta 4')
plt.legend(['Sun','Earth'])
ax.set_xlim([-xmax,xmax])
ax.set_ylim([-ymax,ymax])
ax.set_zlim([-zmax,zmax])
ax.set_xlabel('X[AU]')
ax.set_ylabel('Y[AU]')
ax.set_zlabel('Z[AU]')
plt.show()
"""
