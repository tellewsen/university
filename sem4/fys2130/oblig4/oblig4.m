close all
clear all

%Read file
filnavn = 'tempBlindern10aar.txt';
fileID = fopen(filnavn, 'r');
A = fscanf(fileID, '%d %d %f %f %f',[5,inf]);
minT = A(4,:);
maxT = A(5,:);

n = length(minT);
time = 1:n;

%Plotting
plot(minT,'-r');
hold on;
plot(maxT,'-b');
xlabel('Time[days]');
ylabel('Temperature[C]');
legend('MinT','MaxT')

%Fourier Transform
FMin  = fft(minT);
FMinR = real(FMin);
FMinI = imag(FMin);

FMax  = fft(maxT);
FMaxR = real(FMax);
FMaxI = imag(FMax);

%Plot of Fourier
figure;
plot(time,FMinR,'-r');
hold on
plot(time,FMinI,'-b');
legend('Real','Imaginary')
title('Fourier of MinT')

figure;
plot(time,FMaxR,'-r');
hold on
plot(time,FMaxI,'-b');
legend('Real','Imaginary')
title('Fourier of MaxT')

%Fixing av MaxT
% X:1  , Y:  3.981e+04
% X:11 , Y: -2.099e+04

%Fixing av MinT
%X:1   , Y:  1.321e+04
%X:11  , Y: -1.648e+04
