function bolgeanimering1
clear all;
A = 2;
k = 3;
omega = 8;
N = 1000;
x = linspace(0,20,N);
y = linspace(0,20,N);
p = plot(x,y,'-','EraseMode','xor');
axis([-1 21 -4.5 4.5])
for i=1:200
    t = i*0.01;
    y = A*sin(k*x-omega*t) + A*sin(k*x+omega*t);
    set(p,'XData',x,'YData',y)
    drawnow
    pause(0.02); % For å forsinke fremvisningen
end