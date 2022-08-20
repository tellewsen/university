function Del3
%Konstanter
v       = 340.29; %Speed of sound at sea level [m/s]

%Genererer signalene
N       = 1024*64;
Fs      = 44100;
Fmax    = 1000; % Lowpass frequency max
fsenter = 7000;
fsigma  = 6000;

k01     = HvitStoyGauss(Fs,N,fsenter,fsigma);
k02     = HvitStoyGauss(Fs,N,fsenter,fsigma);
k03     = HvitStoyGauss(Fs,N,fsenter,fsigma);
k04     = HvitStoyGauss(Fs,N,fsenter,fsigma);
k05     = HvitStoyGauss(Fs,N,fsenter,fsigma);
k06     = HvitStoyGauss(Fs,N,fsenter,fsigma);
k07     = HvitStoyGauss(Fs,N,fsenter,fsigma);
k08     = HvitStoyGauss(Fs,N,fsenter,fsigma);
k09     = HvitStoyGauss(Fs,N,fsenter,fsigma);
k10     = HvitStoyGauss(Fs,N,fsenter,fsigma);
k11     = HvitStoyGauss(Fs,N,fsenter,fsigma);

a       = linspace(0,5,1001);
s       = 10;
h       = 400;
dr      = zeros(11,2);
dn      = zeros(11,2);

CCFunk  = zeros(numel(a),1);
%Beregner avstander fra kilder til mottakere

k = 1;
while k <= numel(a) 
    AB = Avstander(a(k),s,h);
    for j=1:2
        for i=1:11
            dr(i,j) = AB(i,j)- h;
            dn(i,j) = round(Fs.*dr(i, j)./v);
        end
    end
    
    %Beregner signalene motatt i A og B
    MottattA = zeros(1,N-max(max(dn)));
    MottattB = zeros(1,N-max(max(dn)));
    for i=1:(N-max(max(dn)))
        MottattA(i) = k01(i+dn(1,1))  + k02(i+dn(2,1)) + k03(i+dn(3,1))...
                    + k04(i+dn(4,1))  + k05(i+dn(5,1)) + k06(i+dn(6,1))...
                    + k07(i+dn(7,1))  + k08(i+dn(8,1)) + k09(i+dn(9,1))...
                    + k10(i+dn(10,1)) + k11(i+dn(11,1));
        
        MottattB(i) = k01(i+dn(1,2))  + k02(i+dn(2,2)) + k03(i+dn(3,2))...
                    + k04(i+dn(4,2))  + k05(i+dn(5,2)) + k06(i+dn(6,2))...
                    + k07(i+dn(7,2))  + k08(i+dn(8,2)) + k09(i+dn(9,2))...
                    + k10(i+dn(10,2)) + k11(i+dn(11,2));
    end
    
    %Lowpass filters
    
    %Mottatt signal i A etter Lowpass
    A2LP = LowPass(MottattA.*MottattA,N,Fs,Fmax);
    %Mottatt signal i B etter Lowpass
    B2LP = LowPass(MottattB.*MottattB,N,Fs,Fmax); 
        
    %Krysskorrelasjon
    CCFunk(k) = CrossCorr(A2LP,B2LP);
    %Counter
    k = k + 1;
end
CCFunk = 1/max(CCFunk).*CCFunk;

snitt = mean(CCFunk((3/4)*numel(a):numel(a)));
maks  = max(CCFunk) - snitt;
linje = maks/2 + snitt;

figure;
plot(a, CCFunk, 'b-', a, linje, 'r-')
title('Plot of CrossCorrelation vs distance between detectors')
xlabel('Distance between detectors [m]')
ylabel('CrossCorrelation[rel.unit]')
end