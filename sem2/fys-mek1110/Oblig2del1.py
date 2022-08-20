from scitools.std import *

#Import measurements from File1
infile = open('motion1.d', 'r')
t = [];  ax = []; ay = []
for line in infile:
    words = line.split()
    t.append(float(words[0]))
    ax.append(float(words[1]))
    ay.append(float(words[2]))
infile.close()

#Turn lists into arrays
t  = array(t)
ax = array(ax)
ay = array(ay)

#Easier notation
n = len(t)
dt = t[2]-t[1]

#Empty arrays for File1 and File2
vx = zeros(n)
vy = zeros(n)
rx = zeros(n)
ry = zeros(n)
for k in range(n-1):
    vx[k+1] = vx[k] + dt*ax[k]
    rx[k+1] = rx[k] + dt*vx[k]

    vy[k+1] = vy[k] + dt*ay[k]
    ry[k+1] = ry[k] + dt*vy[k]

#Plots
figure(1)
plot(rx,ry,legend='Path of object' ,axis=[0,70,0,-35],xlabel='x[m]',ylabel='y[m]',hardcopy='PlotOppgaveE.png')

figure(2)
plot(t, ax, t, ay,legend=['Acc x-axis', 'Acc y-axis'],xlabel='x[m/s**2]',ylabel='y[m/s**2]')

figure(3)
plot(t, vx, t, vy,legend=['Speed x-axis', 'Speed y-axis'],xlabel='x[m/s]',ylabel='y[m/s]')
