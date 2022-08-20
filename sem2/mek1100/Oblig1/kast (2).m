function [ tscaled,xscaled,yscaled ] = kast(x0,y0,v0,theta)
g = 9.81;
t = linspace(0,10,1001);
tscaled = g.*t./(2.*v0.*sin(theta));
xscaled = tscaled;
yscaled = tscaled.*tan(theta) - tscaled.^2.*tan(theta);

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


end

