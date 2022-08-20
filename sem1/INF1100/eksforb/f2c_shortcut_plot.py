from scitools.std import *

F= linspace(-20,120,100)
Capprox = (F-30)/2
Cexact  = (F-32)*5./9


plot(F,Capprox,F,Cexact, legend=['Capprox','Cexact'])

"""
testet funker
"""
