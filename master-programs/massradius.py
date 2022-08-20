import numpy as np
import matplotlib.pyplot as plt
radius,mass = np.loadtxt('Halos/LCDM_vanilla.txt',unpack=True,usecols=(7,8))

plt.plot(radius,mass,'.')
plt.yscale('log')
plt.xlabel(r'r [$Mpc/h$]')
plt.ylabel(r'M [$M_\odot/h$]')
plt.show()

test=mass/(4/3.*np.pi*radius**3)
test=np.log10(test)
plt.plot(radius,test,'.')
plt.show()
