from scitools.StringFunction import StringFunction

f = StringFunction(raw_input('Formula with x: '))

print f(eval(raw_input('x= ')))


"""
Terminal> python user_formula2.py
Formula with x: x+2
x= 1
3.0

Terminal> python user_formula2.py
Formula with x: 3*x+2
x= 1
5
"""

