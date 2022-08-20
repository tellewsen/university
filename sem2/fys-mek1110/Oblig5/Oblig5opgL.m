%Given values
k  = 500;       %N/m
m  = 5;         %kg
h  =.3;         %m
L0 =.5;         %m
my = 0.05;
g  = 9.81;      %m/s^2
G  = -m*g;      %N
%Timesteps
dt= 0.01;       
t = 0:dt:10;    %Timesteps
n= length(t);   %Number of iterations
%Prepare arrays
Fx = zeros(n);  %Vertical Spring force
Fy = zeros(n);  %Horisontal Spring force
N  = zeros(n);  %Normal force
a  = zeros(n);  %Horisontal acceleration
v  = zeros(n);  %Horisontal velocity
x  = zeros(n);  %Horisontal position
y  = zeros(n);  %Vertical position

%Initial values
x(1) = 0.75;    %Initial position

%Euler method
for i = 1:n;
    Fy(i)  = -k*h.*(1-L0./sqrt(x(i).^2+h^2));       %Vertical Spring force
    Fx(i)  = -k.*x(i).*(1-L0./sqrt(x(i).^2+h^2));   %Horisontal Spring force
    N(i)   = -(Fy(i) + G);      %Normal force
    if v(i)== 0;
        Fd(i)= 0;
    else
        Fd(i)  = -(v(i)/abs(v(i)))*my.*N(i);          %Friction force
    end
    %Fd(i)  = 0;
    a(i+1) = (Fx(i)+Fd(i))/m;   %Horisontal acceleration
    v(i+1) = v(i) + dt*a(i+1);    %Horisontal velocity
    x(i+1) = x(i) + dt*v(i+1);    %Horisontal position
end

%Plot
figure(1)
plot(t,x)
xlabel('t[s]')
ylabel('x[m]')
title('Position vs time')

figure(2)
plot(t,v)
xlabel('t[s]')
ylabel('v[m/s]')
title('Velocity vs time')