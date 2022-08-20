#Exercise 6.7 p.332
from scitools.std import *
from sys import *

infile = open('pos.dat', 'r')	#Opens the file (files must be in same folder)
s = float(infile.readline())	#Takes first line and puts it in s
lines = infile.readlines()	#Takes the rest and puts them in lines
infile.close()			#Closes the file

print lines
a = len(lines)			#Used for arrays

x= zeros(a)			#Prepares arrays
y= zeros(a)

for line, i in zip(lines,range(a+1)):	#Fills x and y arrays
    words = line.split()
    x[i] = words[0]
    y[i] = words[1]

vx = zeros(a-1)			#Preparing arrays again
vy = zeros(a-1)
t  = zeros(a-1)

for e in range(a-1):		#Fills t array
    t[e]  = 15*e
    vx[e] = (x[e+1] - x[e]) /s
    vy[e] = (y[e+1] - y[e]) /s

plot(x,y)
wait=raw_input('Showing x versus y, press Return to proceed:')

plot(t,vx,xlabel=t,ylabel=vx)
wait=raw_input('Showing t versus vx,press Return to proceed:')

plot(t,vy,xlabel=t,ylabel=vy)
wait=raw_input('Showing t versus vy,press Return to exit:')

"""
Terminal> python position2velocity.py 
Showing x versus y, press Return to proceed:
Showing t versus vx,press Return to proceed:
Showing t versus vy,press Return to exit:
"""
