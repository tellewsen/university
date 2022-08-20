%Given values
k  = 500;       %N/m
m  = 5;         %kg
h  =.3;         %m
L0 =.5;         %m
my = 0.05;
g  = 9.81;      %m/s^2
G  = -m*g;      %N


x = linspace(0.4,0.75,1000);
Fy  = -k*h.*(1-L0./sqrt(x.^2+h^2));
Fx  = -k.*x.*(1-L0./sqrt(x.^2+h^2));
N   = -(Fy + G);
Fd  = -my.*N;
AllFx = abs(Fx + Fd);
trapz(x,AllFx)

%ans = 25.1074