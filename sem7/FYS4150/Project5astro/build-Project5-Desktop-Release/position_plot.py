import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

"""
This file makes a video of the position of all objects for each time step
"""
#data = np.loadtxt("outputRK4.dat",unpack=True);
data = np.loadtxt("outputverlet.dat",unpack=True);
xdata = data[0]
ydata = data[1]
zdata = data[2]


#Insert number of planets in data
n_planets = 100
t_final = 5
#Insert number of timesteps to create data
timesteps = int(len(xdata)/n_planets)

AU = 1.496E11
ly = 9.4607e15

x = np.zeros((n_planets,timesteps))
y = np.zeros((n_planets,timesteps))
z = np.zeros((n_planets,timesteps))

m = 'o'	#marks objects as dots
c = 'r'	#color of points

for t in range(0,timesteps):
	for i in range(0,n_planets):	
		x[i,t] = xdata[i +t*n_planets]
		y[i,t] = ydata[i +t*n_planets]
		z[i,t] = zdata[i +t*n_planets]


#Plot initial position of all objects
fig = plt.figure(0)
ax = fig.add_subplot(111, projection='3d')
ax.plot(x[:,0],y[:,0],z[:,0],'o')
plt.title('Initial positions of objects')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

xmax = 20.
ymax = 20.
zmax = 20.
#Plot the trajectories of all objects
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
for i in range(0,n_planets):
	ax.plot(x[i,:],y[i,:],z[i,:])
ax.set_xlim([-xmax,xmax])
ax.set_ylim([-ymax,ymax])
ax.set_zlim([-zmax,zmax])
ax.set_xlabel('X[ly]')
ax.set_ylabel('Y[ly]')
ax.set_zlabel('Z[ly]')
ax.set_title('Trajectories of all objects')
plt.show()


fig2 = plt.figure(3)
counter =0 #counter used for saving imagefiles
for t in xrange(0,timesteps,2):
	ax = fig2.add_subplot(111, projection='3d')
	s = ax.scatter(x[:,t], y[:,t], z[:,t],c=c,marker=m)
	s.set_edgecolors = s.set_facecolors = lambda *args:None
	ax.set_xlim([-xmax,xmax])
	ax.set_ylim([-ymax,ymax])
	ax.set_zlim([-zmax,zmax])
	ax.set_xlabel('X[ly]')
	ax.set_ylabel('Y[ly]')
	ax.set_zlabel('Z[ly]')
	ax.set_title(r'Time=%1.3f $\tau_{crunch}$'%((float(t)/timesteps*t_final)))
	#plt.show(False)
	plt.savefig('imagefiles/tmp%05d.png'%counter)
	counter +=1
	fig2.clf()

