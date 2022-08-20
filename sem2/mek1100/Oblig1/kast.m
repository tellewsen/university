function [ output_args ] = kast(x0,y0,v0,theta)
g = 9.81;
t = linspace(0,100);
tstar = g.*t./(2.*v0.*sin(theta));
rx = (v0.*tstar.*cos(theta)-x0).*g/(v0.^2.*sin(2.*theta));
ry = (v0.*tstar.*sin(theta)-0.5.*g.*tstar.^2-y0).*g/(v0.^2.*sin(2.*theta));
plot(rx,ry)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


end

