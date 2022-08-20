close all
% Enkelt eksempelprogram for å vise hvordan fouriertransformasjon
% kan gjennomføres i praksis i Matlab. Eksemplet er en modifikasjon
% av et eksempelprogram på hjelpesidene i Matlab.
Fs = 1000;              % Samplingsfrekvens
delta_t = 1/Fs;         % Tid mellom hver sampling
N = 1024;               % Antall samplinger
t = (0:N-1)*delta_t;    % Tidsvektor
% Lager her et kunstig signal som en sum av et 50 Hz sinussignal
% og en 120 Hz cosinus, pluss legger til et random signal:
x = 0.7*sin(2*pi*50*t) + cos(2*pi*120*t);
x = x + 1.2*randn(size(t));
plot(Fs*t,x)            % Plotting av signalet i tidsbilet
title('Opprinnelig signal (tidsbildet)')
xlabel('tid (millisekunder)')
X = fft(x,N)/N; 

% Fouriertransformasjon
frekv = (Fs/2)*linspace(0,1,N/2); % Frekvensvektor (for plot)

% Plotter bare lengden på frekvenskomponentene i frekvensspekteret.
% Velger å bare ta med frekvenser opp til halve samplingsfrekvensen.
figure;                     % Hindrer overskriving av forrige figur
plot(frekv,2*abs(X(1:N/2))) % Plotter halvparten av fourierspekteret
title('Absolutt-verdier av frekvensspekteret')
xlabel('Frekvens (Hz)')
ylabel('|X(frekv)|')