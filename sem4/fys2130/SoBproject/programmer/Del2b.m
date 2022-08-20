function Del2b

%Genererer signalene k1 og k2
N  = 64*1024; % 2^6*2^10 = 2^16 
Fs = 44100;
%T = N/Fs;
%t = linspace(0,T*(N-1)/N,N);

Fmax = 1000;
fsenter = 10000;
fsigma  = 1000;

k1 = HvitStoyGauss(Fs,N,fsenter,fsigma);
k2 = HvitStoyGauss(Fs,N,fsenter,fsigma);

%Genererer f og g
f = k1+k2;
A2LP = LowPass(f.*f,N,Fs,Fmax);
CrossCorrFunk = zeros(2001,1);
g = zeros(N,1);

for dn=0:2000
    for i=1:60000
        g(i) = k1(i + dn) + k2(i);
    end
    B2LP = LowPass(g.*g,N,Fs,Fmax);
    CrossCorrFunk(dn+1) = CrossCorr(A2LP,B2LP);
end
CrossCorrFunk = 1/max(CrossCorrFunk) .*CrossCorrFunk;

%Regner ut krysskorrelasjon mellom f og g med forskjellige dn
dn = (0:2000);
snitt = mean(CrossCorrFunk((3/4)*2001:2001));
maks = max(CrossCorrFunk) - snitt;
linje = maks/2 + snitt;

%Plotting
figure;
plot(dn,CrossCorrFunk,'b-',dn,linje, 'r-')
title('CrossCorrelation vs dn')
xlabel('dn[m]')
ylabel('CrossCorrelation[rel.unit]')

%Wavelet for forskjellige dn
K    = 48;
fmin = 300;
fmax = 1100;

%f og g1
dn = 0;
g  = zeros(N,1);
for i=1:60000
    g(i) = k1(i + dn) + k2(i);
end
B2LP = LowPass(g.*g,N,Fs,Fmax);
WL1(A2LP,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon
WL1(B2LP,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

%g2
dn = 4;
g  = zeros(N,1);
for i=1:60000
    g(i) = k1(i + dn) + k2(i);
end
B2LP = LowPass(g.*g,N,Fs,Fmax);
WL1(B2LP,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

%g3
dn = 9;
g  = zeros(N,1);
for i=1:60000
    g(i) = k1(i + dn) + k2(i);
end
B2LP = LowPass(g.*g,N,Fs,Fmax);
WL1(B2LP,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

%g4
dn = 26;
g  = zeros(N,1);
for i=1:60000
    g(i) = k1(i + dn) + k2(i);
end
B2LP = LowPass(g.*g,N,Fs,Fmax);
WL1(B2LP,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

end