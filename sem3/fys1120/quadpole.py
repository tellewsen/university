# -*- coding: utf-8 -*-

from pylab import*
close("all")

# Define position and magnitude of electric charges
a = 0.1  #m
q = 1e-9 #C

xPosition = 0.5*array([ a,  -a, -a,   a])
yPosition = 0.5*array([ a,   a, -a,  -a])
charge  = array([q, q, -q, -q])


# Define region size
xMin = -0.4  * a
xMax =  0.4  * a
yMin = -0.4  * a
yMax =  0.4  * a

# Define permittivity of space
epsilon0 = 8.85e-12 


# Initialize arrays to store E field vector components. 
n = 51
Ex = zeros((n,n))
Ey = zeros((n,n))

xPoints = linspace(xMin,xMax,n)
yPoints = linspace(yMin,yMax,n)

# Iterate through each charge on outer loop
# i,j are 'plotting' coordinates
# x,y are 'E field calculation' coordinates
for q in range(4):
  i = 0
  # Now iterate in the region of interest
  for x in xPoints:
    j = 0
    for y in yPoints:
      # Calculate vector components in the charge-to-point direction
      rx = x-xPosition[q]
      ry = y-yPosition[q]

      # Calculate distance r between current point and this charge
      r = sqrt(rx**2+ry**2)

      # Do not count current charge if calculating point on
      # of charge position
      if( r == 0 ):
        continue
      
      
      # calculate unit vector
      rx = rx/r
      ry = ry/r

      # Calculate X and Y contributions for this charge and point, adding
      # it as a vector to the result obtained with previous charges.
      E = (charge[q]/(4*pi*epsilon0*r**2))
      Ex[j,i] = Ex[j,i] + E*rx
      Ey[j,i] = Ey[j,i] + E*ry
      j = j+1 
    i = i+1  


#Calculate and plot E field magnitude
figure(figsize=(8,6))
extent = [xMin,xMax,yMin,yMax]
imshow(sqrt(Ex**2+Ey**2),origin='lower', extent = extent)

#Add colorbar
colorbar()

#Direction of the field
h = quiver(xPoints[::3],yPoints[::3],Ex[::3,::3], Ey[::3,::3]);

# Add title and units, label axes
title("Electric field magnitude [V/m] and direction at each point")
xlabel('X position [m]')
ylabel('Y position [m]')


