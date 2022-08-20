from scitools.std import *

#Constants
c   = 299792458 #m/s
wl0 = 656.3 #m

#Doppler function
def Doppler(wl):
    return c*(wl-wl0)/wl0


#Load data and plot wavelength vs flux for all 10 datafiles
wl,f = loadtxt("spectrum_day0.txt",unpack=True)
subplot(5,2,1)
plot(wl,f)
axis(656.295,656.505,0.6,1.2)

wl,f = loadtxt("spectrum_day67.txt",unpack=True)
subplot(5,2,2)
plot(wl,f)
axis(656.295,656.505,0.6,1.2)

wl,f = loadtxt("spectrum_day133.txt",unpack=True)
subplot(5,2,3)
plot(wl,f)
axis(656.295,656.505,0.6,1.2)

wl,f = loadtxt("spectrum_day200.txt",unpack=True)
subplot(5,2,4)
plot(wl,f)
axis(656.295,656.505,0.6,1.2)

wl,f = loadtxt("spectrum_day267.txt",unpack=True)
subplot(5,2,5)
plot(wl,f)
axis(656.295,656.505,0.6,1.2)

wl,f = loadtxt("spectrum_day333.txt",unpack=True)
subplot(5,2,6)
plot(wl,f)
axis(656.295,656.505,0.6,1.2)

wl,f = loadtxt("spectrum_day400.txt",unpack=True)
subplot(5,2,7)
plot(wl,f)
axis(656.295,656.505,0.6,1.2)

wl,f = loadtxt("spectrum_day467.txt",unpack=True)
subplot(5,2,8)
plot(wl,f)
axis(656.295,656.505,0.6,1.2)

wl,f = loadtxt("spectrum_day533.txt",unpack=True)
subplot(5,2,9)
plot(wl,f)
axis(656.295,656.505,0.6,1.2)

wl,f = loadtxt("spectrum_day600.txt",unpack=True)
subplot(5,2,10)
plot(wl,f)
axis(656.295,656.505,0.6,1.2)
