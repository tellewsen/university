% Energy transport in an Einstein crystal
clear all;
Ly = 40; % System size y-direction
LA = 40; % Length of system A in x direction
LB = 40; % Length of system B in x direction
Lx = LA + LB;
NA = LA*Ly;
NB = LB*Ly;
qA = 3000; % Initial energy in system A
qB = 0; % Initial energy in system B
q = qA + qB; % Total energy - conserved
N = NA + NB; % Total number of oscillators
state = zeros(Lx,Ly); % 2d state matrix
% Generate initial, random states for A and B
for ia = 1:qA
ix = randi(LA,1,1); % Rnd position from 1 to LA
iy = randi(Ly,1,1); % Rnd position from 1 to Ly
state(ix,iy) = state(ix,iy) + 1; % Add energy to this site
end
for ib = 1:qB
ix = randi(LB,1,1)+LA; % Rnd pos from LA+1 to LA+LB
iy = randi(Ly,1,1); % Rnd pos from 1 to Ly
state(ix,iy) = state(ix,iy) + 1; % Add energy to this site
end
% Simulate state development
nstep = 10000000; % nr of simulation steps
EA = zeros(nstep,1); % Energy per oscillator in system A
EB = EA; % Energy per oscillator in system B
for istep = 1:nstep
% Select an oscillator at random
ix1 = randi(Lx,1,1);
iy1 = randi(Ly,1,1);
% Check if this oscillator has non-zero energy
if (state(ix1,iy1)>0)
% Find a random neighbor
dx = 2*randi(2,1,1)-3; % +/-1 with equal prob
ix2 = mod(ix1 + dx-1,Lx)+1; % mod for periodic boundaries
dy = 2*randi(2,1,1)-3; % +/-1 with equal prob
iy2 = mod(iy1 + dy-1,Ly)+1; % mod for periodic boundaries
% Transfer energy from (ix1,iy1) to (ix2,iy2)
state(ix2,iy2) = state(ix2,iy2) + 1;
state(ix1,iy1) = state(ix1,iy1) - 1;
end
if (mod(istep,1000)==0) % Display system at regular intervals
imagesc(state'),colorbar, axis equal, axis tight, drawnow
end
end
