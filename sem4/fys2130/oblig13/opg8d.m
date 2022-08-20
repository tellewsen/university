function opg8d
% Program som foretar en waveletanalyse av et lydsignal som foreligger som
% en wav-fil. Laget for FYS2130 v�ren 2014 av Arnt Inge Vistnes.
% Denne versjonen er fra 13042014.


% PARAMETRE VI M� SETTE: (Parametre for waveletanalysen settes nedenfor!)
N  = 8192;
f1 = 1000; %Hz
f2 = 1600; %Hz
c1 = 1;
c2 = 1.7;
fs = 10000;     %Hz
sigma1 = 0.01;   %s
sigma2 = 0.10;   %s
t1 = 0.15; %s
t2 = 0.5;  %s
T = N/fs; % Total tid lydutsnittet tar (i sek)
t = linspace(0,T*(N-1)/N,N);
h =   c1.*sin(2.*pi.*f1.*t).*exp(-((t-t1)/sigma1).^2)...
    + c2.*cos(2.*pi.*f2.*t).*exp(-((t-t2)/sigma2).^2);

plot(t,h,'-k');
title('Opprinnelig tids-funksjon');
xlabel('Tid (sek)');
ylabel('Signal (rel enhet)');


%************************************************************************
% Foretar s� en effektiv wavelettransformasjon basert p� FFT/IFFT. P� dette
% punktet m� totalt antall punkter N og samplingsfrekvens fs v�re kjent.

% Beregner f�rst FFT av tidsstrengen h
FTsignal = fft(h);

% Plotter frekvensspekteret (absoluttverdier only)
f = linspace(0,fs*(N-1)/N, N);
nmax = floor(N/2); % Plotter bare opp til halve samplingsfrekv.
figure;
plot(f(1:nmax),abs(FTsignal(1:nmax)));
xlabel('Frekvens (Hz)');
ylabel('Relativ intensitet');
title('Vanlig frekvensspektrum');

%************************************************************************
%Ved � stoppes kriptet p� dettestedet, kan vi f� ut FT-spekteret
%og bestemme oss for hvilket frekvensomr�de vi vil bruke for
%waveletanalyse. Skriv s� inn disseverdiene i de neste to linjene
%f�r programmet s� kj�res med WT inne.
%INPUTPARAMETRE vi selv m� velge f�lger her:
fmin=800.0;    %Minimumfrekvensiwaveletanalysen(iHz)
fmax=2000.0;   %Maximumfrekvens
K=100;           %Morlet-wavelet-bredde(kan v�re 6 - 400)

WL1(h,N,fmin,fmax,K,fs);%Kallerp�wavelettransformasjon