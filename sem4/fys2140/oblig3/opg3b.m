%Oblig 3 opg 3b FYS2140 07.02.2014

%Constants
k1   = 0.6;
k2   = 0.7;
hbar = 1;
c    = 1;
A    = 1;
m    = 1;
x    = linspace(-100,100,10000);
w1   = sqrt(k1^2+ 1);
w2   = sqrt(k2^2+ 1);

%Calculations and plot
for t=0:4
    y1 = A*sin(k1.*x-w1.*t);
    y2 = A*sin(k2.*x-w2.*t);
    f = y1+y2;
    plot(x,f)
    xlabel('x')
    ylabel('f(x)')
    title('Plot of f(x) with x=[-100,100] and t=0,1,2,...,30')
    hold on
end