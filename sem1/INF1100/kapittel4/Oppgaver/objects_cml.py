import sys
a= eval(sys.argv[1])
print a, type(a)

"""
Terminal> python objects_cml.py 1
1 <type 'int'>

Terminal> python objects_cml.py 1.1
1.1 <type 'float'>

Terminal> python objects_cml.py "(1,1)"
(1, 1) <type 'tuple'>

Terminal> python objects_cml.py [1,1] 
[1, 1] <type 'list'>

Terminal> python objects_cml.py 2j   
2j <type 'complex'>

Terminal> python objects_cml.py ""This is a string"" 
Traceback (most recent call last):
  File "objects_cml.py", line 2, in <module>
    a= eval(sys.argv[1])
  File "<string>", line 1, in <module>
NameError: name 'This' is not defined

Terminal> python objects_cml.py '"This is a string"'
This is a string <type 'str'>
"""
