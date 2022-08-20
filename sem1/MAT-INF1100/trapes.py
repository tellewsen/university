import time
from math import *

# Funksjonen som skal integreres

def f(x):
	return cos(x)
		
# Integrasjonsfunksjon som bruker trapesmetoden ved aa dele [a,b] i 1, 2, 4, 8, ... 
# Det hele stopper naar vi kommer til 2^M oppdelinger eller estimert relativ feil er mindre enn eps
# Relativ feil estimeres slik som beskrevet i kompendiet

def trapes(f,a,b,eps,M):
	n = 1; h = b-a;
	I = 0.5*h*(f(a)+f(b))
	abserr = abs(I)
	
	j = 1
	while j<M and abserr>=eps*abs(I):
		j = j + 1
		Ip = I

		sum = 0.0; n = 2*n;
		h = h/2.0
		for i in range(1,n):
			x = a + i*h
			sum = sum + f(x)
		
		I = 0.5*h*(f(a)+f(b)+2*sum)
		abserr=abs(I-Ip)
		
	return(I,abserr/I,j)
	
	
# Et kort hovedprogram som setter det hele i gang

a = 0.0; b = 1.0
eps = 1.0e-10; M = 25

[I,err,K] = trapes(f,a,b,eps,M)

# Vi skriver ut resultatene og tidsforbruket.

print "Integral: %1.14f, estimert relativ feil: %e, iterasjoner: %i" %(I, err, K)
print "\nEksakt verdi=%1.14f, eksakt relativ feil=%e" %(sin(b),abs(sin(b)-I)/abs(sin(b)))
