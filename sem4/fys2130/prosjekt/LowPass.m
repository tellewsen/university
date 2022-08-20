function v = LowPass(u, N, Fs, Fmax)
 % Lowpass filter.
 % This program removes all frequencies above a frequency fmax in the
 % frequency spectrum of a signal u. 
 % N is the length of u and Fs is the sampling frequency of u.
 %
 % Returns the real part of the signal after the frequencies have been
 % stripped. It should be noted that the imaginary part of all v are 0
 % before we remove the imaginary part. The only reason we bother to remove
 % it is in hopes of faster computation.
u = fft(u);
nmax =floor(Fmax/(Fs*(N-1))*N*N); 
u(nmax+1:N-nmax+1) = 0;
v = real(ifft(u));
end