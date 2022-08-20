import numpy as np




data = np.loadtxt('Halos/SymA_vanilla.txt',unpack=True)

masses = data[8]


limits = [10**12,10**13.5]
print ((masses>limits[0]) &(masses<limits[1])).sum()
