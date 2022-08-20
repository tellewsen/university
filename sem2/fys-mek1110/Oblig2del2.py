from scitools.std import *

#Import measurements from File1
infile = open('motion1.d', 'r')
t1 = [];  ax1 = []; ay1 = []
for line in infile:
    words = line.split()
    t1.append(float(words[0]))
    ax1.append(float(words[1]))
    ay1.append(float(words[2]))
infile.close()

#Import measurements from File2
infile = open('motion2.d', 'r')
t2 = [];  ax2 = []; ay2 = []
for line in infile:
    words = line.split()
    t2.append(float(words[0]))
    ax2.append(float(words[1]))
    ay2.append(float(words[2]))
infile.close()

#Turn lists into arrays
t1  = array(t1)
ax1 = array(ax1)
ay1 = array(ay1)
t2  = array(t2)
ax2 = array(ax2)
ay2 = array(ay2)

#Easier notation
n1  = len(t1)
dt1 = zeros(len(t1))
n2  = len(t2)
dt2 = zeros(len(t2))

#Empty arrays for File1 and File2
vx1 = zeros(n1)
vy1 = zeros(n1)
rx1 = zeros(n1)
ry1 = zeros(n1)

vx2 = zeros(n2)
vy2 = zeros(n2)
rx2 = zeros(n2)
ry2 = zeros(n2)

#Euler method
for k in range(n1-1):
    dt1[k+1] = t1[k+1] - t1[k]
    vx1[k+1] = vx1[k] + dt1[k]*ax1[k]
    vy1[k+1] = vy1[k] + dt1[k]*ay1[k]
    rx1[k+1] = rx1[k] + dt1[k]*vx1[k]
    ry1[k+1] = ry1[k] + dt1[k]*vy1[k]
for k in range(n2-1):
    dt2[k+1] = t2[k+1] - t2[k]
    vx2[k+1] = vx2[k] + dt2[k]*ax2[k]
    vy2[k+1] = vy2[k] + dt2[k]*ay2[k]
    rx2[k+1] = rx2[k] + dt2[k]*vx2[k]
    ry2[k+1] = ry2[k] + dt2[k]*vy2[k]

#Plots
plot(rx1,ry1, rx2, ry2,axis=[0,70,0,-35], legend=['Path of object(Motion1)','Path of object(Motion2)'],xlabel='x[m]',ylabel='y[m]',hardcopy='PlotOppgaveH.png')


#Exercise G and I
amagnitude1 = zeros(n1)
amagnitude2 = zeros(n2)
for k in range(n1):
    amagnitude1[k] = sqrt(ax1[k]**2+ay1[k]**2)
for k in range(n2):
    amagnitude2[k] = sqrt(ax2[k]**2+ay2[k]**2) 

#Magnitude max searchls
for k in range(n1-1):
    if amagnitude1[k] == max(amagnitude1):
        print 'Magnitude of acceleration max at time %ss with position %sm and value %sm/s**2 for Motion1.d'%(t1[k],rx1[k],amagnitude1[k])
for k in range(n2-1):
    if amagnitude2[k] == max(amagnitude2):
        print 'Magnitude of acceleration max at time %ss with position %sm and value %sm/s**2 for Motion2.d'%(t2[k],rx2[k],amagnitude2[k])
    
