function Del1b
% Lager et tilfeldig signal som ønsket
Fs = 44100;
N = 1024*64;
T = N/Fs;
t = linspace(0,T*(N-1)/N, N);
f = linspace(0,Fs*(N-1)/N,N);
fsenter = 10000;
xx1 = HvitStoyGauss(Fs,N,fsenter,0.5);
xx2 = HvitStoyGauss(Fs,N,fsenter,10);
xx3 = HvitStoyGauss(Fs,N,fsenter,100);
xx4 = HvitStoyGauss(Fs,N,fsenter,1000);
xx5 = HvitStoyGauss(Fs,N,fsenter,5000);

% Ploter frekvens- og tidsbildet
figure;
plot(t,xx1,'-g');
xlabel('Tid (s)');
ylabel('Signal (vilkårlig enhet)');
title('Signal med fsigma=0.5')

figure;
plot(f,abs(fft(xx1)),'-r');
xlabel('Frekvens (Hz)');
ylabel('Frekvenskomponent (rel. verdier)');
title('Fourier for signal med fsigma=0.5')

figure;
plot(t,xx2,'-g');
xlabel('Tid (s)');
ylabel('Signal (vilkårlig enhet)');
title('Signal med fsigma=10')

figure;
plot(f,abs(fft(xx2)),'-r');
xlabel('Frekvens (Hz)');
ylabel('Frekvenskomponent (rel. verdier)');
title('Fourier for signal med fsigma=10')

figure;
plot(t,xx3,'-g');
xlabel('Tid (s)');
ylabel('Signal (vilkårlig enhet)');
title('Signal med fsigma=100')

figure;
plot(f,abs(fft(xx3)),'-r');
xlabel('Frekvens (Hz)');
ylabel('Frekvenskomponent (rel. verdier)');
title('Fourier for signal med fsigma=100')

figure;
plot(t,xx4,'-g');
xlabel('Tid (s)');
ylabel('Signal (vilkårlig enhet)');
title('Signal med fsigma=1000')

figure;
plot(f,abs(fft(xx4)),'-r');
xlabel('Frekvens (Hz)');
ylabel('Frekvenskomponent (rel. verdier)');
title('Fourier for signal med fsigma=1000')

figure;
plot(t,xx5,'-g');
xlabel('Tid (s)');
ylabel('Signal (vilkårlig enhet)');
title('Signal med fsigma=5000')

figure;
plot(f,abs(fft(xx5)),'-r');
xlabel('Frekvens (Hz)');
ylabel('Frekvenskomponent (rel. verdier)');
title('Fourier for signal med fsigma=5000')



% Wavelet analyse
K = 48;             % Morlet-wavelet-bredde (kan være 6 - 400+)

fmin =  9000.0;     % Minimum frekvens i waveletanalysen (i Hz)
fmax = 11000.0;     % Maximum frekvens
WL1(xx1,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

fmin =  9500.0;     % Minimum frekvens i waveletanalysen (i Hz)
fmax = 10500.0;     % Maximum frekvens
WL1(xx2,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

fmin =  9000.0;     % Minimum frekvens i waveletanalysen (i Hz)
fmax = 11000.0;     % Maximum frekvens
WL1(xx3,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

fmin =  8000.0;     % Minimum frekvens i waveletanalysen (i Hz)
fmax = 12000.0;     % Maximum frekvens
WL1(xx4,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

fmin =  4000.0;     % Minimum frekvens i waveletanalysen (i Hz)
fmax = 25000.0;     % Maximum frekvens
WL1(xx5,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

% Autokorrelasjons-analyse
tittel='Autokorrelasjon for Signal1';
AutoCorr(xx1,N,Fs,tittel);

tittel='Autokorrelasjon for Signal2';
AutoCorr(xx2,N,Fs,tittel);

tittel='Autokorrelasjon for Signal3';
AutoCorr(xx3,N,Fs,tittel);

tittel='Autokorrelasjon for Signal4';
AutoCorr(xx4,N,Fs,tittel);

tittel='Autokorrelasjon for Signal5';
AutoCorr(xx5,N,Fs,tittel);
end