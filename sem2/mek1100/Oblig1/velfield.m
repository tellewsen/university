function [X,Y,u,v] = velfield(n)
%setter standardverdi for n til 20
if nargin < 1; 
    n=20;
end

x =linspace(-0.5*pi,0.5*pi,n);
[X,Y] = meshgrid(x,x);

u = cos(X).*sin(Y);
v = -sin(X).*cos(Y);
end