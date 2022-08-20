function u=varmeledning (ustart, u0, u1, N, r)
% løs varmeledningslikninga med startverdi ustart
% og med randkrav u0 til venstre og u1 til høyre,
% utfør N tidsskritt og la r = Delta t / (Delta x)^2 være gitt

u = ustart;
[MM,M]=size(u);
if ((MM > 1) || (M < 3))
  error('Temperaturvektoren har ikke riktig størrelse, må være 1xM hvor M>2')
end
u(1)=u0;
u(M)=u1;
for n = 1:N
  u = r*circshift(u,[0,1]) + (1 - 2*r)*u + r*circshift(u,[0,-1]);
  u(1) = u0;
  u(M) = u1;
end
