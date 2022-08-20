q = [['a','b','c'],['d','e','f'],['g','h']]

print 'Letter a:', q[0][0]
print 'List d,e,f:', q[1]
print 'Last element, h:', q[-1][-1]
print 'Element d:', q[1][0]

"""
q[-1][-2] gives you g since [-1] takes you to the last list, then -2 takes you to the second latest item in that list, which in this case is the first
"""
