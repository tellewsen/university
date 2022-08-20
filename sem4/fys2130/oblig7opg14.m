function bolgeanimasjonX
% Genererer posisjonsarray
delta_x = 0.1;
x = -20:delta_x:20;
n = length(x);
% Genererer posisjoner ved t=0
sigma = 2.0;
u = exp(-(x/(2*sigma)).*(x/(2*sigma))); % Gaussisk form
plot(x,u,’-r’);
figure;
% Genererer div parametre og tidsderivert av utslag ved t=0
v = 0.5;delta_t = 0.1;
faktor = (delta_t*v/delta_x)^2;
dudt = (v/(2*sigma*sigma))*x.*u;
% Angir effektive initialbetingelser:
u_jminus1 = u - delta_t*dudt;
u_j = u;
for t = 1:1000
u_jplus1(2:n-1) = (2*(1-faktor))*u_j(2:n-1) - ...
u_jminus1(2:n-1) + faktor.*(u_j(3:n)+u_j(1:n-2));
% Håndtering av randproblemet, dvs setter u_j(-1) = u_j(n+1)=0
u_jplus1(1) = (2*(1-faktor)).*u_j(1) - u_jminus1(1) + faktor.*u_j(2);
u_jplus1(n) = (2*(1-faktor)).*u_j(n) - u_jminus1(n) + faktor.*u_j(n-1);
plot(u_j);
axis([0 n+1 -0.3 1.2])
drawnow;
u_jminus1 = u_j;
u_j = u_jplus1;
end;