from scitools.std import *

c   = 299792458 #m/s
wl0 = 656.3 #m

def Doppler(wl):
    return c*(wl-wl0)/wl0

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


#By staring at the plot i assume the center of the line
#is at 656.4nm. I insert this into the Doppler formula.


print "v1  = ",Doppler(656.4)
print "v2  = ",Doppler(656.356)
print "v3  = ",Doppler(656.4)
print "v4  = ",Doppler(656.425)
print "v5  = ",Doppler(656.425)
print "v6  = ",Doppler(656.425)
print "v7  = ",Doppler(656.39)
print "v8  = ",Doppler(656.38)
print "v9  = ",Doppler(656.395)
print "v10 = ",Doppler(656.425)
