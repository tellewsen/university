function [] = kast(m,k,u,v,n);
t = linspace(0,n);
g = 9.81;
rx = (m.*u./k.*(1-exp(-k.*t./m)));
ry = -m.*g./k.*t + (m.*v./k + m.^2.*g./k^2).*(1-exp(-k.*t./m));
plot(rx,ry)
xlabel('x')
ylabel('y')
end