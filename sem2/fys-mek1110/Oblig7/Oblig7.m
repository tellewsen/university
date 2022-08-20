% This program also needs the function
% force .m , which must be in the same directory
% as this program
%
% Modify from here -->
N = 3; % nr of balls
m = 0.2; % kg
k = 1000; % N/m
q = 4;
d = 0.05; % m
v0 = 3; % m/s
time = 1; % s
dt = 0.001; %s
% Base variables
n = ceil ( time /dt);
x = zeros (n,N);
v = zeros (n,N);
F = zeros (n,N);
t = zeros (n ,1);
% Initial conditions , equally spaced
for j = 1:N
x(1,j) = d*(j -1);
end
v(1 ,1) = v0;
for i = 1:n -1
% Find force on each block , j
% First , force from block to the left
for j = 2:N
dx = x(i,j) - x(i,j -1);
F(i,j) = F(i,j) + force (dx ,d,k,q);
end
% Second , force from block to the right
for j = 1:N -1
dx = x(i,j +1) - x(i,j);
F(i,j) = F(i,j) - force (dx ,d,k,q);
end
% Euler - Cromer step
%for j = 1:N
%a = F(i,j)/m;
%v(i+1,j) = v(i,j) + a*dt;
%x(i+1,j) = x(i,j) + v(i+1,j)*dt;
%end
%
% The Euler - Cromer step above can also be vectorized
% through the following implementation ( which is faster )
% Euler - Cromer vectorized step
a = F(i ,:) /m;
v(i+1 ,:) = v(i ,:) + a*dt;
x(i+1 ,:) = x(i ,:) + v(i ,:) *dt;
%
t(i+1) = t(i) + dt;
end
% Plot results
for j = 1:N
plot (t,v(:,j));
if j==1
hold on;
end
if j==N
hold off
end
end
xlabel ('t [s]');
ylabel ('v [m/s]');
title('N=3  m=0.2  k=1000  q=4  d=0.05  v0=3  time=1  dt=0.001')

%H) Testet og ser greit ut.
%I) Gir mening ved start, men det som skjer etter det gir ikke mening.
%   v0A =3, v0B,v0C=0. Programmet sier at farten til a blir negativ
%   etterhvert og at farten til b stabiliserer seg ved 0.5m/s.
%   Dette stemmer helt klart ikke med virkeligheten.
%
%J) Plottet blir mer presist. Med dette mener jeg at farten til ball b og c
%   kommer seg naermere null. Det er fortsatt en feil her.
%K) Hvis vi setter variablene slik som jeg har gjort i plottet blir det
%   ganske likt virkeligheten
%
%
%