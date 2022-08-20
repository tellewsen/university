from numpy import exp,linspace
from matplotlib.pyplot import figure,plot,show
T= linspace(0,1E15,1000)
T9 		= T*1E-9		#Convert T to form used in reactions.
T9_star1 	= T9/(1+4.95E-2*T9)	#Other form used.
T9_star2 	= T9/(1+.759*T9)	#Yet another form used.

H_D   	= 	4.01E-15*T9**(-2./3)*exp(-3.380*T9**(-1./3))*(1 +\
                     .123*T9**(1./3) + 1.09*T9**(2./3) + .938*T9)

He3_p	= 	6.04E10*T9**(-2./3)*exp(-12.276*T9**(-1./3))*(1 +\
                     .034*T9**(1./3) - .522*T9**(2./3)-.124*T9 + \
                     .353*T9**(4./3)+ .213*T9**(-5./3))

He3_Be  = 	5.61E6*T9_star1**(5./6)*T9**(-3./2)*exp(-12.826*\
                     T9_star1**(-1./3))

Be_Li   = 	1.34E-10*T9**(-1./2)*(1 - .537*T9**(1./3)  + \
                     3.86*T9**(2./3) +  .0027*T9**-1*exp(2.515E-3*T9**-1))

Li_p	= 	     1.096E9*T9**(-2./3)*exp(-8.472*T9**(-1./3)) - \
                     4.830E8*T9_star2**(5./6)*T9**(-3./2)*exp(-8.472*\
                     T9_star2**(-1./3)) + 1.06E10*T9**(-3./2)\
                     *exp(-30.442*T9**-1)

esum = H_D + He3_p + He3_Be + Be_Li + Li_p
figure()
plot(T,esum)
show()
