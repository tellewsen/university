function y = HvitStoyGauss(Fs,N,fsenter,fsigma)

%y = zeros(N,1);
%T = N/Fs;
%t = linspace(0,T*(N-1)/N,N);
f = linspace(0,Fs*(N-1)/N, N);
%nsenter = floor(N*fsenter/(Fs*(N-1)/N));
%nsigma = floor(N*fsigma/(Fs*(N-1)/N));
gauss = exp(-(f-fsenter).*(f-fsenter)/(fsigma*fsigma));
%figure;
%plot(f,gauss,'-k'); % For sjekking!

ampl = rand(N,1);
ampl = ampl.*transpose(gauss);
%figure;
%plot(f,ampl,'-g'); % For sjekk

faser = rand(N,1);
faser = faser*2*pi;
y = ampl.*(cos(faser) + 1i.*sin(faser));

%Fikser spekteret slik at kriterier for folding er ok
y(1) = 0 + 0i;
for i=1:N/2-1
    y(N/2+1+i) = conj(y(N/2+1-i));
end

%Plotter fourier spekteret til hvit stoy
% figure;
% plot(f,y)
% title('Fourier spekteret til hvit stoy')

y = ifft(y);

% %test
% figure;
% plot(f,real(y))
% title('Real del av signalet')
% figure;
% plot(f,imag(y))
% title('Imaginaer del av signalet')

y = real(y);
end