from scitools.std import *

#Load data for day0
wl,F = loadtxt("spectrum_day0.txt",unpack=True)

#Define functions
def Fmodel(wl,Fmin,sigma,lambdacenter):	
    return Fmax + (Fmin - Fmax)*exp(-(wl-lambdacenter)**2 / (2*sigma**2))

def Delta(wl, Fmin, sigma, lambdacenter):
    return sum((F - Fmodel(wl,Fmin,sigma,lambdacenter))**2)

#Arrays and constants
n 	     = 20				#Number of variables for each interval
Fmax  	     = 1.				#Given in assignment
Fmin         = linspace(  0.7 , 0.9   ,n)	#Interval1
sigma        = linspace(  0.015, 0.025  ,n)	#Interval2
lambdacenter = linspace(656.37, 656.43,n)	#Interval3

#Best fit model
best_Fmin         = Fmin[0]			#Save first Fmin as best
best_sigma        = sigma[0]			#Save first sigma as best
best_lambdacenter = lambdacenter[0]		#Save first lamdacenter as best
best_delta	  = 50 				#Some number that is larger than the first delta

#For loop to calculate best combination of variables
for i in range(n):
    for j in range(n):
        for k in range(n):
            delta = Delta(wl, Fmin[i], sigma[j], lambdacenter[k])
	    if delta <= best_delta:
	        best_Fmin = Fmin[i]
		best_sigma = sigma[j]
		best_lambdacenter = lambdacenter[k]
		best_delta = delta
#Plotting
plot(wl,F)							#Plot Wavelength vs Flux
hold("on")							#Make plots appear in same window
plot(wl,Fmodel(wl,best_Fmin,best_sigma,best_lambdacenter))	#Plot Wavelength vs Fluxmodel
axis(656.295,656.505,0.65,1.2)					#Fix axis
xlabel("Wavelength[m]"), ylabel("Flux[W/m^2]")			#Name axes
savefig("Day0.png")						#Save plot to file

#Prints the fits to terminal
print "Fmin:    ",best_Fmin
print "Sigma:   ",best_sigma
print "Lcenter: ",best_lambdacenter
print "Delta:   ",best_delta

"""
Terminal:
thorae@vanir ~/privat/sem3/ast1100/oblig3 $ python Andreas_Ellewsen_Oblig3_del2.py 
Fmin:     0.805263157895
Sigma:    0.0202631578947
Lcenter:  656.398421053
Delta:    2.46356272361
"""
