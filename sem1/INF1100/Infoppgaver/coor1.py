h = 0.01

x_list = []

for i in range (0, 101):
    x = 1 + i*h
    x_list.append(x)
for e in x_list:
    print "%.2f" %(e)
