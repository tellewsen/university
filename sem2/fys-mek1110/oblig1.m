T=20;
dt=0.01;
n=ceil(T/dt);
t=zeros(n,1);
Fc=zeros(n,1);
Fv=zeros(n,1);
D=zeros(n,1);
v=zeros(n,1);
a=zeros(n,1);
x=zeros(n,1);
F=zeros(n,1);
rho=1.293;
Cd=1.2;
A=0.45;
w=0;
tc=0.67;
fc=488;
fv=25.8;
i=1;
m=80;

while (i<n-1)&&(x(i)<100)
D(i)=(1/2)*A*(1-0.25*exp(-(t(i)./tc).^2))*rho*Cd*(v(i).^2);
Fc(i)=fc*exp(-(t(i)./tc).^2);
Fv(i)=fv.*v(i);
F(i)=400;
a(i)=(F(i)+Fc(i)-Fv(i)-D(i))/m;
v(i+1)=v(i) + dt*a(i);
x(i+1) = x(i) + dt*v(i);
t(i+1) = t(i) + dt;
i = i + 1;
end
subplot(2,2,1)
plot(t(1:i-1),F(1:i-1))
xlabel('t [s]')
ylabel('F [N]')
subplot(2,2,2)
plot(t(1:i-1),Fc(1:i-1))
xlabel('t [s]')
ylabel('Fc [N]')
subplot(2,2,3)
plot(t(1:i-1),Fv(1:i-1))
xlabel('t [s]')
ylabel('Fv [N]')
subplot(2,2,4)
plot(t(1:i-1),D(1:i-1))
xlabel('t [s]')
ylabel('D [N]')
display(t(i))