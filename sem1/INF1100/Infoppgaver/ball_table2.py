g  = 9.81
v0 = 10
T  = 2*v0/g
n  = 81
dt = T/(n-1)

t = [i*dt for i in range(0, n)]
y = [v0*t_i - 0.5*g*t_i**2 for t_i in t]

for t_i, y_i in zip(t, y):
    print '%.1f  %.4f' % (t_i, y_i)



