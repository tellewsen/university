#Exercise 2.2 p. 79

F = 0                                                #Startvalue for F(Farenheit degrees)
dF= 10                                               #Value used in last line of while loop
while F <= 100:                                      #Checks that F has not passed 100 and stops when it has
    C = (F-32.0)*(5.0/9.0)                           #Calculates C (Celsius Degrees)
    Ca= (F-30.0)/2.0                                 #Calculates Ca (C approximate)
    print "F= %3.0f  C= %6.2f  Ca = %6.2f" %(F,C,Ca) #Prints each of the lines in the table
    F = F + dF                                       #Increases F by dF (In this case by 10 each run)


"""
Terminal> python f2c_approx_table.py 
F=   0  C= -17.78  Ca = -15.00
F=  10  C= -12.22  Ca = -10.00
F=  20  C=  -6.67  Ca =  -5.00
F=  30  C=  -1.11  Ca =   0.00
F=  40  C=   4.44  Ca =   5.00
F=  50  C=  10.00  Ca =  10.00
F=  60  C=  15.56  Ca =  15.00
F=  70  C=  21.11  Ca =  20.00
F=  80  C=  26.67  Ca =  25.00
F=  90  C=  32.22  Ca =  30.00
F= 100  C=  37.78  Ca =  35.00
"""
