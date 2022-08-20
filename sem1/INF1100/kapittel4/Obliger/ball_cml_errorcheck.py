import sys 				#Needed for sys.argv
g = 9.81 				#Gravity..

try:					#Tries to get v0 from commandline
    v0= float(sys.argv[1])		
except IndexError:			#Goes here if you don't give v0 on commandline
    v0= float(raw_input('v0= '))	#Asks for v0 and turns answer into float
except ValueError:			#Goes here if you write something else than a number
    print 'v0 must be a number'		#Tells you what you did wrong
    v0= float(raw_input('v0= '))	#Asks for v0 and turns answer into float

try:
    t = float(sys.argv[2])		#Same thing as above for t
except IndexError:
    t = float(raw_input('t= '))
except ValueError:
    print 't must be a number'
    t= float(raw_input('t= '))

if 0< t <2*v0/g:			#Tests if t is betwen 0 and 2*v0/g and exits the program if it is
    y = v0*t - 0.5*g*t**2		#The equation
    print 'y = %g' %y			#Printstatement
else:
    print 'The t you input is not between 0 and 2*v0/g and this program does not allow this.. sorry bro.. exiting...'
    sys.exit()



"""
Terminal> Obliger $ python ball_cml_errorcheck.py 100 1
y = 95.095

Terminal> python ball_cml_errorcheck.py 100 2
y = 180.38

Terminal> python ball_cml_errorcheck.py 100 3
y = 255.855

Terminal> python ball_cml_errorcheck.py 1 50 
The t you input is not between 0 and 2*v0/g and this program does not allow this.. sorry bro.. exiting...

Terminal> python ball_cml_errorcheck.py
v0= 1 
t= 50
The t you input is not between 0 and 2*v0/g and this program does not allow this.. sorry bro.. exiting...

Terminal> python ball_cml_errorcheck.py
v0= 100
t= 1
y = 95.095
"""
