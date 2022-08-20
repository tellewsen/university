%Given values
k  = 500;       %N/m
m  = 5;         %kg
h  =.3;         %m
L0 =.5;         %m
my = 0.05;
g  = 9.81;      %m/s^2
G  = -m*g;      %N

%Timesteps
dt = 0.001;      
t  = 0:dt:10;    %Timesteps
n  = length(t);  %Number of iterations

%Prepare arrays
Fx    = zeros(n,1);  %Vertical Spring force
Fy    = zeros(n,1);  %Horisontal Spring force
N     = zeros(n,1);  %Normal force
Kin   = zeros(n,1);  %Kinetic energy
Fnetx = zeros(n,1);  %Horisontal sum of forces
Pot   = zeros(n,1);  %Potential energy
a     = zeros(n,1);  %Horisontal acceleration
v     = zeros(n,1);  %Horisontal velocity
x     = zeros(n,1);  %Horisontal position

%Initial values
x(1)  = 0.75;        %Initial position

%Euler method
for i = 1:n;
    Fy(i)  = -k*h.*(1-L0./sqrt(x(i).^2+h^2));     %Vertical Spring force
    Fx(i)  = -k.*x(i).*(1-L0./sqrt(x(i).^2+h^2)); %Horisontal Spring force
    N(i)   = -(Fy(i) + G);                        %Normal force
    if v(i) == 0;
        Fd(i)= 0;
    else
        Fd(i)  = -(v(i)/abs(v(i)))*my.*N(i);      %Friction force
    end
    %Fd(i)  = 0;
    a(i+1) = (Fx(i)+Fd(i))/m;   %Horisontal acceleration
    v(i+1) = v(i) + dt*a(i+1);  %Horisontal velocity
    x(i+1) = x(i) + dt*v(i+1);  %Horisontal position
    Kin(i+1) = 1/2*m.*v(i+1).^2;%Kinetic energy
    Fnetx(i+1) = Fx(i) + Fd(i); %Horisontal sum of forces
end

Pot = 25.1074 + cumtrapz(x,-Fnetx);

figure(1)
plot(x,Kin,'r')
title('Kinetic Energy vs position')
axis([-.8,.8,0,25])
xlabel('x[m]')
ylabel('Kinetic Energy[J]')
figure(2)
plot(x,Pot,'b')
axis([-.8,.8,0,25])
title('Potential Energy vs position')
xlabel('x[m]')
ylabel('Potential Energy[J]')