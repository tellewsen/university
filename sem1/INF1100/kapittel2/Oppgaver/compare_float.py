a = 1/947.0*947
b = 1
if abs(a-b) < 1.0e-16:
    print 'Wrong result!'

"""
Terminal> python compare_float.py 
"""

"Instead of checking a != b we check if abs(1-b) is less than a really small number, and decide that they really are the same"
