v0 = 1. #m/s
g= 9.81 #m/s**2
step = (2.*v0)/(g*10.)

t_list = []
y_list = []

print '  t          y'
for i in range (0,11):
    t = i*step
    t_list.append(t)
    y = v0*t - 0.5*g*t**2
    y_list.append(y)
    print '%.3f %10.3f' %(t, y)
