function Del2a


%Genererer signalene k1 og k2
N  = 64*1024; % 2^6*2^10 = 2^16 
Fs = 44100;
fsenter = 10000;
fsigma  = 1000;

k1 = HvitStoyGauss(Fs,N,fsenter,fsigma);
k2 = HvitStoyGauss(Fs,N,fsenter,fsigma);

%Genererer f og g for alle forskyvninger dn og
%regner ut krysskorrelasjon mellom f og g for forskjellige dn
f = k1+k2;
g = zeros(N,1);
CrossCorrFunk = zeros(2001, 1);

for dn=0:2000
    for i=1:60000
        g(i) = k1(i + dn) + k2(i);
    end
    CrossCorrFunk(dn+1) = 1/N.*CrossCorr(f,g);
end

dn = (0:2000);
snitt = mean(CrossCorrFunk( (3/4)*2001:2001) );
maks = max(CrossCorrFunk) - snitt;
linje = maks/2 + snitt;


figure;
plot(dn,CrossCorrFunk,'b-',dn,linje, 'r-')
title('CrossCorrelation vs dn')
xlabel('dn[rel.unit]')
ylabel('CrossCorrelation[rel.unit]')


% Wavelet analyse
K = 48;         % Morlet-wavelet-bredde (kan være 6 - 400+)
fmin =  8000.0;  % Minimum frekvens i waveletanalysen (i Hz)
fmax = 12000.0; % Maximum frekvens
WL1(f,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

%g1
dn = 0;
g = zeros(N,1);
for i=1:60000
    g(i) = k1(i + dn) + k2(i);
end
WL1(g,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

%g2
dn = 4;
g = zeros(N,1);
for i=1:60000
    g(i) = k1(i + dn) + k2(i);
end
WL1(g,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

%g3
dn = 9;
g = zeros(N,1);
for i=1:60000
    g(i) = k1(i + dn) + k2(i);
end
WL1(g,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

%g4
dn = 13;
g = zeros(N,1);
for i=1:60000
    g(i) = k1(i + dn) + k2(i);
end
WL1(g,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

%g5
dn = 82;
g = zeros(N,1);
for i=1:60000
    g(i) = k1(i + dn) + k2(i);
end
WL1(g,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

end