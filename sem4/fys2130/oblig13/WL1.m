function WL1(h,N,fmin,fmax,K,fs)
% Funksjon for kontinuerlig diskret Morlet waveletanalyse av en array h.
% Array har N elementer samplet ved samplingsfrekvensen fs.
% Morlet waveletanalysen skal gjennomføres i frekvensintervallet
% [fmin,fmax] med K-parameteren K. Se detaljer i et kapittel i læreboka
% Svingninger og bølger, Arnt Inge Vistnes, 2014.

FTsignal = fft(h);
T = N/fs; % Total tid lydutsnittet tar (i sek)
t = linspace(0,T*(N-1)/N,N);
f = linspace(0,fs*(N-1)/N, N);

% Beregner # analysefrekvenser, skriver til skjerm, klargjør frekvensene
M = floor(log(fmax/fmin) / log(1+(1/(8*K)))) + 1;
AntallFrekvenserIAnalyse = M
ftrinn = (fmax/fmin)^(1/(M-1));
f_analyse = fmin;

% Allokerer plass til waveletdiagrammet og array for lagring av frekvenser
WLdiagram = zeros(M,N);
fbrukt = zeros(1,M);

%Løkke over alle frekvenser som inngår i analysen
for jj=1:M
    faktor=(K/f_analyse)*(K/f_analyse);
    FTwl=exp(-faktor*(f-f_analyse).*(f-f_analyse));
    FTwl=FTwl-exp(-K*K)*exp(-faktor*(f.*f));%Litekorreksjonsledd
    FTwl=2.0*FTwl;      %Faktor(ulikevalg!)
    %Beregnersåenhellinjeiwaveletdiagrammetiénjafs!
    %WLdiagram(jj,:)=abs(ifft(FTwl.*transpose(FTsignal))); %Ettalternativ
    %WLdiagram(jj,:)=sqrt(abs(ifft(FTwl.*(FTsignal))));%Ettannet
    WLdiagram(jj,:)=sqrt(abs(ifft(FTwl.*transpose(FTsignal))));%Et til
    %HUSK AT DU KANSKJE MÅ SETTE PÅ TRASPOSE igjen her
    
    
    %Brukerdensistevariantenforåfåsvakepartierbedresynlig
    fbrukt(jj)=f_analyse; %Lagrerfrekvensenesomfaktiskerbrukt
    f_analyse=f_analyse*ftrinn; %Beregnernestefrekvens
end;
%Redusererfilstørrelsevedåfjernemyeavoverflødiginformasjonitid.
%Dettegjøreskunforatfilstørrelsenpåplotteneskalblihåndterbar.
P = floor((K*fs)/(24*fmax));%Tallet24kanendresvedbehov
TarBareMedHvertXITid = P
NP = floor(N/P);
AntallPktITid = NP
for jj=1:M
    for ii=1:NP
        WLdiagram2(jj,ii)=WLdiagram(jj,ii*P);
        tP(ii)=t(ii*P);
    end;
end;
%Foretaenmarkeringiplottetforåviseområdermedrandproblemer
maxverdi=max(WLdiagram2);
mxv=max(maxverdi);
for jj=1:M
    m=floor(K*fs/(P*pi*fbrukt(jj)));
    WLdiagram2(jj,m)=mxv/2;
    WLdiagram2(jj,NP-m)=mxv/2;
end;

%Plotterwaveletdiagrammet
figure;
imagesc(tP,log10(fbrukt),WLdiagram2,'YData',[1 size(WLdiagram2,1)]);
set(gca,'YDir','normal');
xlabel('Tid(sek)');
ylabel('Log10(frekvensiHz)');
%title(’WaveletPowerSpektrum’);%Velgdennenårdeteraktuelt,
title('Sqrt(WaveletPowerSpektrum)');%mendennenårsqrtblirbrukt
colorbar('location','southoutside');
end