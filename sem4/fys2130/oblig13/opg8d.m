function opg8d
% Program som foretar en waveletanalyse av et lydsignal som foreligger som
% en wav-fil. Laget for FYS2130 våren 2014 av Arnt Inge Vistnes.
% Denne versjonen er fra 13042014.


% PARAMETRE VI MÅ SETTE: (Parametre for waveletanalysen settes nedenfor!)
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
% Foretar så en effektiv wavelettransformasjon basert på FFT/IFFT. På dette
% punktet må totalt antall punkter N og samplingsfrekvens fs være kjent.

% Beregner først FFT av tidsstrengen h
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
%Ved å stoppes kriptet på dettestedet, kan vi få ut FT-spekteret
%og bestemme oss for hvilket frekvensområde vi vil bruke for
%waveletanalyse. Skriv så inn disseverdiene i de neste to linjene
%før programmet så kjøres med WT inne.
%INPUTPARAMETRE vi selv må velge følger her:
fmin=800.0;    %Minimumfrekvensiwaveletanalysen(iHz)
fmax=2000.0;   %Maximumfrekvens
K=100;           %Morlet-wavelet-bredde(kan være 6 - 400)

WL1(h,N,fmin,fmax,K,fs);%Kallerpåwavelettransformasjon