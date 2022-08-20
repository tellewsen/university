for i in range(2, 4):
    print i
    for j in range(i-1, i+1):
        for k in range(j-1, j):
            if i != j:
                print j, k

"""
i = 2
print 2 
j = 1
k = 0 
	i != j 
print 1 0
j = 2
	i /= j 
k ...

i = 3 
print 3 
j = 2 
k = 1
print 2 1 
j = 3
k ...


ender altså opp med:

2
1 0
3
2 1
"""
