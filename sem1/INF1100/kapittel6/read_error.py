#Exercise 6.11 p.333

#To get this program to work you need to run lnsum.py and store the output as lnsum.dat

import sys
from scitools.std import *

def lnread(filename):
    infile = open(filename)
    lines = infile.readlines()[24:]
    a = len(lines)
    epsilon = zeros(a)
    error = zeros(a)
    n = zeros(a)
    for line, e in zip(lines,range(a)):
        epsilon[e] = int(line[9])*10**int(line[11:14])
        error[e] = float(line[29:33])*10**float(line[34:37])
        n[e] = line.split('=')[1]
    return {'epsilon':epsilon, 'error':error, 'n':n}


lnsum = lnread('lnsum.dat')

plot(lnsum['n'],lnsum['error'],log='y')
hold('on')
plot(lnsum['n'],lnsum['epsilon'],log='y')
legend('n vs error','n vs epsilon')


"""
Terminal> python read_error.py
"""
