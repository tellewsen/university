M = 10;
N = 20;
ustart = zeros(1,M);

for t = [0,0.1,0.2,0.3,0.4,0.5]
    dt = 0.1/N;
    x  = linspace(0,1,M);
    dx = 1/(M-1);
    r  = dt/dx.^2;
    u = varmeledning(ustart,0,1,N,r);
    plot(x,u)
    hold on
    ustart = u;
end