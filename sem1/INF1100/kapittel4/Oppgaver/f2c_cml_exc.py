import sys
try:
   F = float(sys.argv[1])
except:
   F = float(raw_input('F= '))
C = (F-32.0)*(5.0/9.0)
print '%g degrees Fahrenheit = %g degress Celsius' %(F,C)

"""
Terminal> python f2c_cml_exc.py 20
20 degrees Fahrenheit = -6.66667 degress Celsius

Terminal> python f2c_cml_exc.py   
F= 20
20 degrees Fahrenheit = -6.66667 degress Celsius

Terminal> python f2c_cml_exc.py 1
1 degrees Fahrenheit = -17.2222 degress Celsius

Terminal> python f2c_cml_exc.py  
F= 1
1 degrees Fahrenheit = -17.2222 degress Celsius
"""
