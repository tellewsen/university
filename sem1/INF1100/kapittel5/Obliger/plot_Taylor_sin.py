#Exercise 5.21 p. 
from scitools.std import *	#imports everything needed(pi,factorial,sin,plot and so on)

n= 70			#used later
x= pi*2.0		#used later
S= 0			#needs to exist for next stuff to work

for j in range(n):	#first part of exercise
    S = S + ((-1.)**(j)*x**(2.*j+1)) / factorial(2.*j+1)


def S(x,n):		#Defines the function for approximating sin with a sum function
    S=0
    for j in range(n):
        S = S + ((-1.)**(j)*x**(2.*j+1)) / factorial(2.*j+1)
    return S

x=linspace(0,4*pi,201)	#array from 0 -> 4pi with 200 intervals


plot(x,sin(x))		#plots sin of those 200 intervals
axis([0,4*pi,-1,1])	#locks the axis so we can see what's going on
hold('on')		#makes the plots appear in the same window

for n in (1,2,3,6,12):
	plot(x,S(x,n))	#plots the approximation of sin for n=1,2,3,6,12

legend(['sin(x)','S(x,1)','S(x,2)','S(x,3)','S(x,6)','S(x,12)'])	#Marks the lines in the plot
title('Approximations of Sin(x)')	#Gives you a title on top of the plot
xlabel('x')				#Labels x axis as x
ylabel('S(x)')				#Labels y axis as S(x)

"""
Terminal>
*plot appears*
"""
