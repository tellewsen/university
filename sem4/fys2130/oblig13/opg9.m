function opg9
% Program som foretar en waveletanalyse av et lydsignal som foreligger som
% en wav-fil. Laget for FYS2130 våren 2014 av Arnt Inge Vistnes.
% Denne versjonen er fra 13042014.

% Navn på noen lydfiler (her angitt sammen med egnede parametre for analyse)
 c = 'Svarttrost2.wav'; % Nstart 17000, 64 k, 1500-8000 Hz, K12-96
% c = 'bokfink.wav';    % Nstart 34 000, 128 k, K 48, 2000-8000 Hz (minst)
% c = 'gjok.wav';        % Nstart 12000, 450-750 Hz + noe støy ved > 4000 Hz
% c = 'tilElise.wav';    % Nstart 100, 128 k, K 48, 90 - 800 Hz

% PARAMETRE VI MÅ SETTE: (Parametre for waveletanalysen settes nedenfor!)
N = 1024*128;   % Lengde på data-utsnitt (lydfil) (helst 2^n) Her: "64k"
nstart = 17000; % Startpunkt for utsnitt fra lydfil

% Leser passe utsnitt av lydfil, spiller den av og plotter tidsbildet
nslutt = nstart+N-1;
[y, fs] = wavread(c, [nstart nslutt]); % Les array y(N,2) fra fil
% 'fs' er vanligvis 44100 (samplingsfrekvens ved CD kvalitet)
h = zeros(N,1); % Plukker ut bare én kanal fra stereosignalet lest
h = y(:,1);
sound(h,fs); % Spiller av utsnittet som er brukt
T = N/fs; % Total tid lydutsnittet tar (i sek)
t = linspace(0,T*(N-1)/N,N);
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
fmin=1500.0;    %Minimumfrekvensiwaveletanalysen(iHz)
fmax=8000.0;   %Maximumfrekvens
K=12;           %Morlet-wavelet-bredde(kan være 6 - 400)

WL1(h,N,fmin,fmax,K,fs);%Kallerpåwavelettransformasjon