clear all
% Define position and magnitude of electric charges
a = 0.1;  %m
q = 1e-9; %C

xPosition = 0.5*[ a,  -a, -a,   a];
yPosition = 0.5*[ a,   a, -a,  -a];
charge  = [q, q, -q, -q];


% Define region size
xMin = -0.4  * a;
xMax =  0.4  * a;
yMin = -0.4  * a;
yMax =  0.4  * a;

% Define permittivity of space
epsilon0 = 8.85e-12;


% Initialize arrays to store E field vector components. 
n = 100;
Ex = zeros(n,n);
Ey = zeros(n,n);


xPoints = linspace(xMin,xMax,n);
yPoints = linspace(yMin,yMax,n);


% Iterate through each charge on outer loop
% i,j are 'plotting' coordinates
% x,y are 'E field calculation' coordinates
for q = 1:4
  i = 0;
 
  % Now iterate in the region of interest
  for x = xPoints
    i = i+1;
    j = 0;
    for y = yPoints
      j = j+1;
    
      % Calculate vector components in the charge-to-point direction
      rx = x-xPosition(q);
      ry = y-yPosition(q);

      % Calculate distance r between current point and this charge
      r = sqrt(rx^2+ry^2);

      % Do not count current charge if calculating point on
      % of charge position
      if( r == 0 )
        continue;
      end

      % Calculate unit vector
      rx = rx/r;
      ry = ry/r;

      % Calculate X and Y contributions for this charge and point, adding
      % it as a vector to the result obtained with previous charges.
      E = (charge(q)/(4*pi*epsilon0*r^2));

      Ex(j,i) = Ex(j,i) + E*rx;
      Ey(j,i) = Ey(j,i) + E*ry;

    end
  end
end


% Calculate and plot E field magnitude
imagesc(xPoints,yPoints, sqrt(Ex.^2+Ey.^2));

% Add colorbar
colorbar

% Flip Y-axis of image (by default this 
% type of plot is upside-down)
set(gca,'YDir','normal');

% Hold plot to overlap quiver on top
hold on 

% Overlap E field vectors, scaling arrows by a factor
% of 2 so the direction can be seen clearly, draw arrow
%for every fifth point
h = quiver(xPoints(1:5:end),yPoints(1:5:end), Ex(1:5:end,1:5:end), Ey(1:5:end,1:5:end),3);
set(h, 'LineWidth',1,'Color' , [0 0 0]);
hold off

% Add title and units, label axes
title('Electric field magnitude [V/m] and direction at each point', 'fontsize', 14);
xlabel('X position [m]', 'fontsize', 14);
ylabel('Y position [m]', 'fontsize', 14);


