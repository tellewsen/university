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

y = v0*t - 0.5*g*t**2			#The equation
print y					#Prints answer to terminal

"""
Terminal> python ball_cml_qa.py one two
v0 must be a number
v0= 1
t must be a number
t= 2
-17.62

Terminal> python ball_cml_qa.py 1 2    
-17.62

Terminal> python ball_cml_qa.py    
v0= 1
t= 2
-17.62
"""
