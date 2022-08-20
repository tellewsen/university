#Exercise 5.12 p.246
from scitools.std import *
import numpy as np

Flist = np.array(range(-20,121))#Makes array with points between -20 and 120.)
Celist = np.zeros(len(Flist))	#Prepares empty array for later
Cilist = np.zeros(len(Flist))	#Same as above

#Define functions
def Cexact(f):
    return (f-32.0)*5.0/9

def Cinexact(f):
    return (f-30.0)/2.0

#Fill lists with functionvalues
for i in range(len(Flist)):
    Celist[i] = Cexact(i)

for i in range(len(Flist)):
    Cilist[i] = Cinexact(i)

plot(Flist,Celist,'r-') 	#Plots Celist, red line
hold('on')			#Keeps first plot
plot(Flist,Cilist,'b-') 	#Plots Cilist, blue line
xlabel('Fahrenheit')	
ylabel('Celsius')
legend('exact','inexact')	#shows which line is which in the plot

"""
Terminal>  python f2c_shortcut_plot.py 
*plot appears*
"""
