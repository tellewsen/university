import random,sys

N = int(sys.argv[1])
r = []
head= 0

for i in range(N):
    r.append(random.random())
    if r[i] <= 0.5:
        print 'head'
    else:
        print 'tail'

for i in range(len(r)):
   if r[i] <= 0.5:
      head +=1

print 'Drew %g times, got head %g times' %(N,head)
