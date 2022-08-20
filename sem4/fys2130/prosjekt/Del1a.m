function Del1a
% Lager et tilfeldig signal som ønsket
Fs = 44100;
N = 1024*64;
T = N/Fs;
t = linspace(0,T*(N-1)/N, N);
f = linspace(0,Fs*(N-1)/N,N);
fsigma = 50;
xx1 = HvitStoyGauss(Fs,N,100,fsigma);
xx2 = HvitStoyGauss(Fs,N,1000,fsigma);
xx3 = HvitStoyGauss(Fs,N,10000,fsigma);

% Velger å plotte tids- og frekvensbildet her
figure;
plot(t,xx1,'-g');
xlabel('Tid (s)');
ylabel('Signal (vilkårlig enhet)');
title('Signal med fsenter=100')

figure;
plot(f,abs(fft(xx1)),'-r');
xlabel('Frekvens (Hz)');
ylabel('Frekvenskomponent (rel. verdier)');
title('Fourier for signal med fsenter=100')

figure;
plot(t,xx2,'-g');
xlabel('Tid (s)');
ylabel('Signal (vilkårlig enhet)');
title('Signal med fsenter=1000')

figure;
plot(f,abs(fft(xx2)),'-r');
xlabel('Frekvens (Hz)');
ylabel('Frekvenskomponent (rel. verdier)');
title('Fourier for signal med fsenter=1000')

figure;
plot(t,xx3,'-g');
xlabel('Tid (s)');
ylabel('Signal (vilkårlig enhet)');
title('Signal med fsenter=10000')

figure;
plot(f,abs(fft(xx3)),'-r');
xlabel('Frekvens (Hz)');
ylabel('Frekvenskomponent (rel. verdier)');
title('Fourier for signal med fsenter=10000')

% Wavelet-analyse
K = 48;         % Morlet-wavelet-bredde (kan være 6 - 400+)

fmin = 30.0;    % Minimum frekvens i waveletanalysen (i Hz)
fmax = 170.0;   % Maximum frekvens
WL1(xx1,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

fmin = 800.0;   % Minimum frekvens i waveletanalysen (i Hz)
fmax = 1200.0;  % Maximum frekvens
WL1(xx2,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

fmin = 9000.0;  % Minimum frekvens i waveletanalysen (i Hz)
fmax = 11000.0; % Maximum frekvens
WL1(xx3,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

%Autokorrelasjons-analyse
tittel='Autokorrelasjon for Signal1';
AutoCorr(xx1,N,Fs,tittel);

tittel='Autokorrelasjon for Signal2';
AutoCorr(xx2,N,Fs,tittel);

tittel='Autokorrelasjon for Signal3';
AutoCorr(xx3,N,Fs,tittel);
end
