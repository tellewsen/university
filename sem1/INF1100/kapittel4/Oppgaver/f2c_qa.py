F= float(raw_input('F: '))
C = (F-32.0)*(5.0/9.0)
print '%g degrees Fahrenheit = %g degress Celsius' %(F,C)


"""
Terminal> python f2c_qa.py 
F: 24.5
24.5 degrees Fahrenheit = -4.16667 degress Celsius

Terminal> python f2c_qa.py 
F: 37
37 degrees Fahrenheit = 2.77778 degress Celsius

Terminal> python f2c_qa.py 
F: 36
36 degrees Fahrenheit = 2.22222 degress Celsius

Terminal> python f2c_qa.py 
F: 35
35 degrees Fahrenheit = 1.66667 degress Celsius

Terminal> python f2c_qa.py 
F: 34
34 degrees Fahrenheit = 1.11111 degress Celsius

Terminal> python f2c_qa.py 
F: 33
33 degrees Fahrenheit = 0.555556 degress Celsius

Terminal> python f2c_qa.py 
F: 32
32 degrees Fahrenheit = 0 degress Celsius
"""
