from scitools.all import *

#Constants
M             = 1
n             = 1000
dtau          = 0.01
Spin          = 38*M
Energypermass = 8.03

#Arrays
r      = zeros(n)
phi    = zeros(n)
x      = zeros(n)
y      = zeros(n)

#Initial Conditions
r[0]      = 20
phi[0]    = 0

# Calculate path of spaceship
for i in range(n-1):
   if r[i] >= 2:
      dphi = Spin/r[i]**2*dtau 
      dr   = -sqrt( (Energypermass**2) - (1+ ((Spin/r[i])**2)*(1-2*M/r[i])))*dtau
      phi[i+1] = phi[i] + dphi
      r[i+1] = r[i] + dr
   else:
      print "Final phi: ", phi[i]/pi*180
      break

#Convert to Cartesian:
for i in range(n):
   x[i] = r[i]*cos(phi[i])
   y[i] = r[i]*sin(phi[i])

#Plot path of spaceship
plot(x,y,'ro')
hold('on')

#Plot eventhorizon of blackhole
rs    = zeros(100)
xs    = zeros(len(rs))
ys    = zeros(len(rs))
rs[:] = 2
phis = linspace(0,2*pi,len(rs))
for i in range(len(rs)):
    xs[i] = rs[i]*cos(phis[i])
    ys[i] = rs[i]*sin(phis[i])
plot(xs,ys,'b-')
xlabel("x-axis[Units of M]")
ylabel("y-axis[Units of M]")
hardcopy("oblig11.png")
