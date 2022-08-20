function TestDel1
% Lager et tilfeldig signal som ønsket
Fs = 44100;
N = 1024*64;
T = N/Fs;
t = linspace(0,T*(N-1)/N, N);
f = linspace(0,Fs*(N-1)/N,N);
fsenter = 10000;
fsigma = 2500;
xx = HvitStoyGauss(Fs,N,fsenter,fsigma);


% % Velger å plotte tids- og frekvensbildet her
% figure;
% plot(t,xx,'-g');
% xlabel('Tid (s)');
% ylabel('Signal (vilkårlig enhet)');
% 
% figure;
% plot(f,abs(fft(xx)),'-r');
% xlabel('Frekvens (Hz)');
% ylabel('Frekvenskomponent (rel. verdier)');
% 
% % Wavelet analyse
% fmin = 7000.0;  % Minimum frekvens i waveletanalysen (i Hz)
% fmax = 13000.0; % Maximum frekvens
% K = 48;         % Morlet-wavelet-bredde (kan være 6 - 400+)
% WL1(xx,N,fmin,fmax,K,Fs); % Kaller på wavelettransformasjon

% Autokorrelasjons-analyse
%AutoCorr(xx,N,Fs);
end
