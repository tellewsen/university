%Given values
k  = 500; %m
m  = 5;   %N/m
h  =.3;   %m
L0 =.5;   %m

%Timesteps
dt= 0.01;       
t = 0:dt:10;
n= length(t);        %Number of iterations

%Prepare arrays
Fx = zeros(n);
a = zeros(n); %Acceleration in horisontal direction
v = zeros(n); %Velocity in horisontal direction
x = zeros(n); %Position in horisontal direction
y = zeros(n);

%Initial values
x(1) = 0.65;    %Initial position
y(1) = 0.3;     %Initial position

%Euler method
for i = 1:n;
Fx(i) = -k.*x(i).*(1-L0./sqrt(x(i).^2+h^2)); %Spring force horisontal
y(i) = h;           
a(i+1) = Fx(i)/m;   %Acceleration in horisontal direction
v(i+1) = v(i) + dt*a(i+1); %Velocity in horisontal direction
x(i+1) = x(i) + dt*v(i+1); %Position in horisontal direction
end

%Plot
figure(1)
plot(t,x)
xlabel('t[s]')
ylabel('x[m]')
title('Position x at time t')

hold on
figure(2)
plot(t,v)
xlabel('t[s]')
ylabel('v[m/s]')
title('Velocity at time t')
