import numpy as np
from scipy import special

# constants
keV = 8.61734e-5  # Boltzmann constant in eV/deg
k_erg = 1.380658e-16  # Boltzmann constant in erg K
h = 6.62607e-27  # Planck constant in erg s
c = 2.99792e10  # Speed of light in cm/s
m_e = 9.109390e-28  # electron mass in grams


# SSA2
# definition of functions
def part_func_e(temp):
    # Partition function U_r
    k_ev_t = keV * temp
    chi_ion = np.array([7, 16, 31, 51])
    u = np.zeros(4)
    for r in range(4):
        for s in range(chi_ion[r]):
            u[r] = u[r] + np.exp(-s / k_ev_t)
    return u  # returns all the values of u array


def boltzmann_e(temp, r, s):
    # Boltzmann distribution n_r,s/N_r
    u = part_func_e(temp)
    rel_nrs = 1. / u[r - 1] * np.exp(-(s - 1) / (keV * temp))
    return rel_nrs


def saha_e(temp, el_press, ion_stage):
    k_ev_t = keV * temp
    k_erg_t = k_erg * temp
    el_dens = el_press / k_erg_t
    chi_ion = np.array([7, 16, 31, 51])
    u = part_func_e(temp)
    u = np.append(u, 2)  # append element to array
    saha_const = (2 * np.pi * m_e * k_erg_t / (h ** 2)) ** (3. / 2) * 2 / el_dens
    n_stage = np.zeros(5)
    n_stage[0] = 1.
    for r in range(4):
        n_stage[r + 1] = n_stage[r] * saha_const * u[r + 1] / u[r] * np.exp(-chi_ion[r] / k_ev_t)
    n_total = np.sum(n_stage)
    n_stage_rel = n_stage / n_total
    return n_stage_rel[ion_stage - 1]


def saha_bolt_e(temp, el_press, ion, level):
    return saha_e(temp, el_press, ion) * boltzmann_e(temp, ion, level)


def saha_bolt_h(temp, el_press, level):
    k_ev_t = keV * temp
    k_erg_t = k_erg * temp
    el_dens = el_press / k_erg_t
    # energy levels and weights for hydrogen
    nr_levels = 100  # Pick a reasonable cutoff value for partition function
    g = np.zeros((2, nr_levels))  # declaration weights
    chi_exc = np.zeros((2, nr_levels))  # declaration excitation energies

    for s in range(nr_levels):
        g[0, s] = 2. * (s + 1.) ** 2  # statistical weights
        chi_exc[0, s] = 13.598 * (1. - 1. / (s + 1.) ** 2.)  # excitation weights
    g[1, 0] = 1.
    chi_exc[1, 0] = 0.

    # Partition functions
    u = np.zeros([2])
    for s in range(nr_levels):
        u[0] = u[0] + g[0, s] * np.exp(-chi_exc[0, s] / k_ev_t)
    u[1] = g[1, 0]

    # Saha
    saha_const = (2. * np.pi * m_e * k_erg_t / (h * h)) ** 1.5 * 2 / el_dens
    n_stage = np.zeros(2)
    n_stage[0] = 1.
    n_stage[1] = n_stage[0] * saha_const * u[1] / u[0] * np.exp(-13.598 / k_ev_t)
    n_total = np.sum(n_stage)  # sum both stages = total hydrogen density

    # Boltzmann
    n_level = n_stage[0] * g[0, level - 1] / u[0] * np.exp(-chi_exc[0, level - 1] / k_ev_t)
    n_level_rel = n_level / n_total  # Fraction of total hydrogen density

    # Test block
    # for s in range(6):
    #     print s+1, g[0, s], chi_exc[0, s], g[0, s]*np.exp(-chi_exc[0, s]/k_ev_t)
    # for s in range(0, nr_levels,10):
    #     print s+1, g[0, s], chi_exc[0, s], g[0, s]*np.exp(-chi_exc[0, s]/k_ev_t)

    return n_level_rel


def part_func_ca(temp):
    # Partition function U_r
    k_ev_t = keV * temp
    chi_ion = np.array([6.113, 11.871, 50.91, 67.15])
    u = np.zeros(4)
    for r in range(4):
        for s in range(int(chi_ion[r])):
            u[r] = u[r] + np.exp(-s / k_ev_t)
    return u  # returns all the values of u array


def boltzmann_ca(temp, r, s):
    # Boltzmann distribution n_r,s/N_r
    k_ev_t = keV * temp
    u = part_func_ca(temp)
    rel_nrs = 1. / u[r - 1] * np.exp(-(s - 1) / k_ev_t)
    return rel_nrs


def saha_ca(temp, el_press, ion_stage):
    k_ev_t = keV * temp
    k_erg_t = k_erg * temp
    el_dens = el_press / k_erg_t
    chi_ion = np.array([6.113, 11.871, 50.91, 67.15])
    u = part_func_ca(temp)
    u = np.append(u, 2)  # append element to array
    saha_const = (2. * np.pi * m_e * k_erg_t / (h * h)) ** (3. / 2) * 2. / el_dens
    n_stage = np.zeros(5)
    n_stage[0] = 1.
    for r in range(4):
        n_stage[r + 1] = n_stage[r] * saha_const * u[r + 1] / u[r] * np.exp(-chi_ion[r] / k_ev_t)
    n_total = np.sum(n_stage)
    n_stage_rel = n_stage / n_total
    return n_stage_rel[ion_stage - 1]


def saha_bolt_ca(temp, el_press, ion_stage, level):
    return saha_ca(temp, el_press, ion_stage) * boltzmann_ca(temp, ion_stage, level)


# SSA3
# Planck function wavelength
def planck(temp, wavelength):
    b_planck = 2.*h*c*c/(wavelength**5) / (np.exp((h*c) / (wavelength * k_erg * temp)) - 1)
    return b_planck


# Voigt function
def voigt(a, u):
    # Calculates the voigt function for values u and a
    z = u + 1.j*a
    return special.wofz(z).real


# Profile function
def profile(a, tau0, u, t_s, wav, t_l):
    # t_s = 5700
    # t_l = 4200
    # wav = 5000e-8
    intensity = np.zeros(u.size)
    for i in range(u.size):
        tau = tau0 * voigt(a, u[i])
        intensity[i] = planck(t_s, wav) * np.exp(-tau) + planck(t_l, wav) * (1. - np.exp(-tau))
    return intensity
