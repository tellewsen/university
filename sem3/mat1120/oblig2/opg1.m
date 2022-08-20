B = [1 -5;2 3];
[x,y] = meshgrid(-15:1:15);

F1 = x-5.*y;
F2 = 2.*x+3.*y;
quiver(x,y,F1,F2)
C1 = [.1,.2,.3,.4,.5];
C2 = [.2,.3,.4,.5,.6];
t = linspace (0,1,100);
for j=1:5
y1 = exp(2.*t).*(-1./2.*C1(j).*cos(3.*t)-3./2.*C1(j).*sin(3.*t)...
    +C2(j)./2.*sin(-3.*t) - 3./2.*C2(j)*cos(-3.*t));
y2 = exp(2.*t).*(C1(j).*cos(3.*t) + C2(j).*sin(-3.*t));
hold on
plot(y1,y2)
end

%For å finne den generelle løsningen må vi radredusere A:
A = [(-3*i-1), -5; 2,(-3*i +1) ];
rref(A)

%%Gir oss [1, 1/2 - 3/2*i; 0, 0]