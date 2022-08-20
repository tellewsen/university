import sys
F= float(sys.argv[1])
C = (F-32.0)*(5.0/9.0)
print '%g degrees Fahrenheit = %g degress Celsius' %(F,C)

"""
Terminal> python f2c_cml.py 5
5 degrees Fahrenheit = -15 degress Celsius

Terminal> python f2c_cml.py 500
500 degrees Fahrenheit = 260 degress Celsius

Terminal> python f2c_cml.py 0  
0 degrees Fahrenheit = -17.7778 degress Celsius

Terminal> python f2c_cml.py 32
32 degrees Fahrenheit = 0 degress Celsius
"""
