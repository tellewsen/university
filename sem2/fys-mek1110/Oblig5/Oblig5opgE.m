k  = 500; %m
m  = 5;   %N/m
h  =.3;   %m
L0 =.5;   %m

x = linspace(-.75,.75,100);
Fx = -k.*x.*(1-L0./sqrt(x.^2+h^2));

plot(x,Fx)
title('Springforce vs position')
xlabel('x[m]')
ylabel('Horisontal springforce[N]')