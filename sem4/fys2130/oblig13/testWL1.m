function testWL1
% Program som foretar en waveletanalyse av et lydsignal som foreligger som
% en wav-fil. Laget for FYS2130 v�ren 2014 av Arnt Inge Vistnes.
% Denne versjonen er fra 13042014.

% Navn p� noen lydfiler (her angitt sammen med egnede parametre for analyse)
% c = 'Svarttrost2.wav'; % Nstart 17000, 64 k, 1500-8000 Hz, K12-96
  c = 'bokfink.wav';    % Nstart 34 000, 128 k, K 48, 2000-8000 Hz (minst)
% c = 'gjok.wav';        % Nstart 12000, 450-750 Hz + noe st�y ved > 4000 Hz
% c = 'tilElise.wav';    % Nstart 100, 128 k, K 48, 90 - 800 Hz

% PARAMETRE VI M� SETTE: (Parametre for waveletanalysen settes nedenfor!)
N = 1024*128;   % Lengde p� data-utsnitt (lydfil) (helst 2^n) Her: "64k"
nstart = 34000; % Startpunkt for utsnitt fra lydfil

% Leser passe utsnitt av lydfil, spiller den av og plotter tidsbildet
nslutt = nstart+N-1;
[y, fs] = wavread(c, [nstart nslutt]); % Les array y(N,2) fra fil
% 'fs' er vanligvis 44100 (samplingsfrekvens ved CD kvalitet)
h = zeros(N,1); % Plukker ut bare �n kanal fra stereosignalet lest
h = y(:,1);
sound(h,fs); % Spiller av utsnittet som er brukt
T = N/fs; % Total tid lydutsnittet tar (i sek)
t = linspace(0,T*(N-1)/N,N);
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
fmin=1500.0;    %Minimumfrekvensiwaveletanalysen(iHz)
fmax=10000.0;   %Maximumfrekvens
K=6;           %Morlet-wavelet-bredde(kan v�re 6 - 400)

WL1(h,N,fmin,fmax,K,fs);%Kallerp�wavelettransformasjon