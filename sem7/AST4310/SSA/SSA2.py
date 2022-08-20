# importing useful libraries (you may need more)
from functions import *
import matplotlib.pyplot as plt  # plotting package
from matplotlib import rc
rc('font', **{'family': 'serif'})  # This is for Latex writing

# Part2
# Compute partition function for T=5000,10000,20000
for temp in [5000, 10000, 20000]:
    print part_func_e(temp)

# Compute distribution for T=5000,
distribution_5k = np.zeros(11)
distribution_10k = np.zeros(11)
distribution_20k = np.zeros(11)
for s in range(1, 11):
    distribution_5k[s - 1] = boltzmann_e(5000, 1, s)
    distribution_10k[s - 1] = boltzmann_e(10000, 1, s)
    distribution_20k[s - 1] = boltzmann_e(20000, 1, s)
print distribution_5k
print distribution_10k
print distribution_20k

# Compute saha_e for different temp
for r in range(1, 6):
    print saha_e(20000, 1e3, r)
print "----"
for r in range(1, 6):
    print saha_e(10000, 1e3, r)

# Compute Saha_bolt_E for different temperatures
for s in range(1, 6):
    print saha_bolt_e(5000, 1e3, 1, s)
print "----"
for s in range(1, 6):
    print saha_bolt_e(20000, 1e3, 1, s)
print "----"
for s in range(1, 6):
    print saha_bolt_e(10000, 1e3, 2, s)
print "----"
for s in range(1, 6):
    print saha_bolt_e(20000, 1e3, 4, s)

# Plotting the population vs temperature for s=1 and P_e = 131
temp = np.arange(0, 30001, 1000)
pop = np.zeros((5, 31))
for T in np.arange(1, 31):
    for r in np.arange(1, 5):
        pop[r, T] = saha_bolt_e(temp[T], 131., r, 1)
label_list = ['ground stage', 'first ion stage', 'seconds ion stage', 'third ion stage']

# Ground state plot
plt.figure(0)
for i in range(1, 5):
    plt.plot(temp, pop[i, :], label=label_list[i - 1])

plt.xlabel('temperature', size=14)
plt.ylabel('population', size=14)
plt.yscale('log')
plt.ylim([1e-3, 1.1])
plt.legend(loc='best')
plt.title('Population vs Temperature for s=1')


# Plotting the population vs temperature for s=2
temp = np.arange(0, 30001, 1000)
pop = np.zeros((5, 31))
for T in np.arange(1, 31):
    for r in np.arange(1, 5):
        pop[r, T] = saha_bolt_e(temp[T], 131., r, 2)
label_list = ['ground stage', 'first ion stage', 'seconds ion stage', 'third ion stage']

# Ground state plot
plt.figure(1)
for i in range(1, 5):
    plt.plot(temp, pop[i, :], label=label_list[i - 1])

plt.xlabel('temperature', size=14)
plt.ylabel('population', size=14)
plt.yscale('log')
plt.ylim([1e-3, 1.1])
plt.legend(loc='best')
plt.title('Population vs Temperature for s=2')

# Plotting the population vs temperature for s=4
temp = np.arange(0, 30001, 1000)
pop = np.zeros((5, 31))
for T in np.arange(1, 31):
    for r in np.arange(1, 5):
        pop[r, T] = saha_bolt_e(temp[T], 131., r, 4)
label_list = ['ground stage', 'first ion stage', 'seconds ion stage', 'third ion stage']

# ground state plot
plt.figure(2)
for i in range(1, 5):
    plt.plot(temp, pop[i, :], label=label_list[i - 1])

plt.xlabel('temperature', size=14)
plt.ylabel('population', size=14)
plt.yscale('log')
plt.ylim([1e-3, 1.1])
plt.legend(loc='best')
plt.title('Population vs Temperature for s=4')
plt.show()

# Printing hydrogen levels
saha_bolt_h(5000, 1e2, 1)  # this works as it should

# Solar Ca+K versus Ha:line strength
temp = np.arange(1000, 20001, 100)
CaH = np.zeros(temp.shape)
Ca_abundance = 2e-6
for i in range(0, len(temp)):
    NCa = saha_bolt_ca(temp[i], 1e2, 2, 1)
    NH = saha_bolt_h(temp[i], 1e2, 2)
    CaH[i] = NCa * Ca_abundance / NH

plt.figure(3)
plt.plot(temp, CaH, label=r'strength ratio Ca$^+$K /H$\alpha$')
plt.yscale('log')
plt.xlabel(r'temperature $T / K$', size=14)
plt.ylabel(r'Ca II K / H$\alpha$', size=14)
plt.legend(fontsize=14)
plt.title('Ca/H Ratio versus temperature')

print 'Ca/H ratio at 5000 K = ', CaH[np.argwhere(temp == 5000)][0][0]

# solar Ca+K versus Ha: temperature sensitivity
temp = np.arange(2000, 12001, 100)
dNCadT = np.zeros(temp.shape)
dNHdT = np.zeros(temp.shape)
dT = 1.
for i in range(101):
    NCa = saha_bolt_ca(temp[i], 1e2, 2, 1)
    NCa2 = saha_bolt_ca(temp[i] - dT, 1e2, 2, 1)
    dNCadT[i] = (NCa - NCa2) / (dT * NCa)
    NH = saha_bolt_h(temp[i], 1e2, 2)
    NH2 = saha_bolt_h(temp[i] - dT, 1e2, 2)
    dNHdT[i] = (NH - NH2) / (dT * NH)

plt.figure(4)
plt.plot(temp, np.abs(dNHdT), label=r'H')
plt.plot(temp, np.abs(dNCadT), ls='--', label=r'Ca$^+$K')
plt.yscale('log')
plt.title('Relative population changes')
plt.ylim(1e-9, 1e1)

NCa = np.zeros(temp.shape)
NH = np.zeros(temp.shape)
for i in range(101):
    NCa[i] = saha_bolt_ca(temp[i], 1e2, 2, 1)
    NH[i] = saha_bolt_h(temp[i], 1e2, 2)

plt.plot(temp, NH / np.amax(NH), label='rel. pop H')
plt.plot(temp, NCa / np.amax(NCa), ls='--', label=r'rel. pop Ca$^+$K')
plt.xlabel(r'temperature $T/K$', size=14)
plt.ylabel(r'$\left| \left( \Delta n(r,s) / \Delta T\right) /  n(r,s) \right|$', size=20)
plt.legend(loc=4, fontsize=12)


# Hot stars vs cold stars
"Find at which temperature the hydrogen in stellar photo spheres with P_e = 100 is about 50% ionized."
for T in np.arange(7000, 10001, 100):
    print T, saha_bolt_h(T, 1e2, 1)
temp = np.arange(6000, 13001, 1e2)
nH = np.zeros(temp.shape)
for i in range(len(temp)):
    nH[i] = saha_bolt_h(temp[i], 1e2, 1)

plt.figure(6)
plt.plot(temp, nH)
plt.xlabel(r'temperature $T/K$', size=14)
plt.ylabel(r'neutral hydrogen fraction', size=14)
plt.title('Fraction of neutral hydrogen in stellar photo spheres')
plt.show()
