import matplotlib.pyplot as plt
import numpy as np
N = 1000
k = 1.38E-23#Boltzmans constant
T_divided_by_theta = np.linspace(0.001, 10, N)
Z_sum = np.zeros(N)
def Partition_function(x):
	Z = 1.
	for j in range(1, 1000):
		Z+= (2*j+1)*np.exp(-j*(j+1)/x)
	return Z
for i in range(N):
	Z_sum[i] = Partition_function(T_divided_by_theta[i])

f = 1 + 3*np.exp(-2/T_divided_by_theta)

plt.figure()
plt.plot(T_divided_by_theta, Z_sum, '-r',
         T_divided_by_theta, f, '-b', 
         T_divided_by_theta, T_divided_by_theta, '-g')
plt.xlabel('T/Theta')
plt.ylabel('Rotational partition function Z')
plt.legend(['Z', 'T << Theta', 'T >> Theta'])
plt.show()


print('Difference in ', Z_sum[-1], T_divided_by_theta[-1])