print '-------------------'
F = 0
dF= 10
while F <= 100:
    C = (F-32.0)*(5.0/9.0)
    print "%.0f, %.2f" %(F,C)
    F = F + dF
print '-------------------'