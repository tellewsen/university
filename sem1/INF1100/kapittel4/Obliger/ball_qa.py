v0 = float(raw_input('v0= ')); t = float(raw_input('t= ')); g = 9.81
y = v0*t - 0.5*g*t**2
print y

"""
Terminal> python ball_qa.py 
v0= 0
t= 5
-122.625

Terminal> python ball_qa.py 
v0= 0
t= 0
0.0

Terminal> python ball_qa.py 
v0= 50
t= 0
0.0

Terminal> python ball_qa.py 
v0= 1
t= 1
-3.905

The only thing I've done here is replace the values for v0 and t with raw_input, which is converted to float afterwards"
"""
