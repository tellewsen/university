n=3
for i in range(-1, n):
    if i != 0:
        print i
for i in range(1, 13, 2*n):
    for j in range(n):
        print i, j
for i in range(1, n+1):
    for j in range(i):
        if j:
            print i, j
for i in range(1, 13, 2*n):
    for j in range(0, i, 2):
            for k in range(2, j, 1):
                b=i>j>k
                if b:
                    print i, j, k






"""
i=-1
print -1
i = 0
i = 1
print 1
i= 2
print 2 
i = 3
print 3
i= 1

skip 2. forloop

i=2
print 2,1
i=3
print 3,1
print 3,2

loop 4:
i=1
j=0
k=2
b= 1>0>2

i=7
j=0
k=2
b= 7>0>2
j=2
k=2
b=7>2>2
j=4
k=2
b=7>4>2
print 7,4,2
osv...



Output:
-1
1
2
3
1 0
1 1
1 2
7 0
7 1
7 2
2 1
3 1
3 2
7 4 2
7 4 3
7 6 2
7 6 3
7 6 4
7 6 5
"""
