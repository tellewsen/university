q = [['a','b','c'],['d','e','f'],['g','h']]

print "a     =", q[0][0]
print "list2 =", q[1]
print "h     =", q[-1][-1]
print "d     =", q[1][0]
print ""
print "-----------fun part below----------"
print ""
print "this should make e:   ", q[-2][-2]
print "a can also be written:", q[-3][-3]
print "b in a weird way:     ", q[0][-2]
print "b the regular way:    ", q[0][1]
print ""
print "-----------eve more fun below------"
print ""

for i in q:
    for j in range(len(i)):
         print i[j]
