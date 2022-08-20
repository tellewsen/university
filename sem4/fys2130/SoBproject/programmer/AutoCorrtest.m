function AutoCorrtest
N= 2^16;
Fs= 44100;

fid = fopen('testsignal.bin','r');
yy = fread(fid, N, 'double');
fclose(fid);

AutoCorr(yy,N,Fs);
end