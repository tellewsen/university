function F = AutoCorr(xx,N,Fs,tittel)
M = N/2;
%Lager funksjonen
F = zeros(N,1);

for j=1:M
    F(j) = sum( xx(1:M).*xx(j:M+j-1) );
    
end
F = 1/F(1) .*F;
T = N/Fs;
tau = linspace(0,T*(N-1)/N, N);

%Plotting
F = abs(F);
Fmax = max(F);

figure;
plot(tau(1:N/2),F(1:N/2),'b-'...
    ,tau(1:N/2),Fmax/2  ,'r-') %Funksjonen

xlabel('Tid[s]')
ylabel('Autokorrelasjonsfunk[rel.enh]')
title(tittel)
hold off
