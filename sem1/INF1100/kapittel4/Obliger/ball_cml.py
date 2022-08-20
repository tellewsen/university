import sys #needed for sys.argv
v0 = float(sys.argv[1]); t = float(sys.argv[2]); g = 9.81
y = v0*t - 0.5*g*t**2
print y


"""
Terminal> python ball_cml.py 1 3
-41.145

Terminal>  python ball_cml.py 0 5
-122.625

Terminal>  python ball_cml.py 0 0
0.0

Terminal>  python ball_cml.py 0 1
-4.905

Terminal>  python ball_cml.py 1 0
0.0

Terminal>  python ball_cml.py 1 50
-12212.5

Same thing as ball_qa.py 
Changed raw_input in v0 and t to sys.argv[1] and sys.argv[2]
"""
