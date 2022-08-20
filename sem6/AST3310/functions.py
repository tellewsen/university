from math import exp
import numpy as np
import sys

# Constants needed
m_u = 1.66053892173E-27  # [kg] Atomic mass unit
sigma = 5.67E-8  # [W*m**-2*K**-4] Stefan-Boltzmann's constant
k = 1.3806488E-23  # [m**2*kg*s**-2*K**-1] Boltzmann's constant
c = 2.99792458E8  # [m*s**-1] Speed of lights in vacuum
a = 4. * sigma / c
N_A = 6.0221413E23  # [] Avogadro's number
MeVtoJ = 1.60217657E-13  # [J] Conversion factor from MeV to Joules

# Energy output of reactions
Q_pp = (.15 + 1.02) * MeVtoJ  # [J]
Q_dp = 5.49 * MeVtoJ  # [J]
Q_33 = 12.86 * MeVtoJ  # [J]
Q_34 = 1.59 * MeVtoJ  # [J]
Q_7e = .05 * MeVtoJ  # [J]
Q_71prime = 17.35 * MeVtoJ  # [J]

# Defining functions for pressure,temperature and density
"""This part solves the equation for Pressure, assuming 
the equation of state to be that of an ideal gas"""


def pressure(rho, t, mu):
    p_rad = a / .3 * t ** 4  # Radiation pressure
    p_g = rho * k * t / (mu * m_u)  # Gas pressure, assuming ideal gas
    p = p_g + p_rad  # Total Pressure
    return p


def density(p, t, mu):
    rho = (p * mu * m_u) / (k * t)  # Ideal gas
    return rho


# Computes number of atoms of the different elements.
def atom_numbers(rho, x, y3, y4, z_73li, z_74be):
    n_e = rho * (1 + x) / (2. * m_u)  # Number of electrons
    n_d = 0  # Number of deuterium atoms
    n_p = rho * x / m_u  # Number of Hydrogen atoms
    n_3 = rho * y3 / (3. * m_u)  # Number of Helium 3 atoms
    n_4 = rho * y4 / (4. * m_u)  # Number of Helium 4 atoms
    n_li = rho * z_73li / (7. * m_u)  # Number of Beryllium 7 atoms
    n_be = rho * z_74be / (7. * m_u)  # Number of Lithium 7 atoms
    numbers = [n_e, n_d, n_p, n_3, n_4, n_li, n_be]
    return numbers


# Reactions
"""In this part the reaction rates of the processes in the PPI and PPII \
chain are calculated. I neglect the other processes since I choose to \
look at a star with low temperature at it's core, and in that case those \
processes are very slow compared to PPI and PPII."""


def reactions(t):
    t9 = t * 1E-9  # Convert T to form used in reactions.
    t9_star1 = t9 / (1 + 4.95E-2 * t9)  # Other form used.
    t9_star2 = t9 / (1 + .759 * t9)  # Yet another form used.

    h_d = 4.01E-15 * t9 ** (-2. / 3) * exp(-3.380 * t9 ** (-1. / 3)) * \
        (1 + .123 * t9 ** (1. / 3) + 1.09 * t9 ** (2. / 3) + .938 * t9)

    he3_p = 6.04E10 * t9 ** (-2. / 3) * exp(-12.276 * t9 ** (-1. / 3)) * (1 + .034 * t9 ** (1. / 3) - .522 * t9 ** (
            2. / 3) - .124 * t9 + .353 * t9 ** (4. / 3) + .213 * t9 ** (-5. / 3))

    he3_be = 5.61E6 * t9_star1 ** (5. / 6) * t9 ** (-3. / 2) * exp(-12.826 * t9_star1 ** (-1. / 3))

    be_li = 1.34E-10 * t9 ** (-1. / 2) * (1 - .537 * t9 ** (1. / 3) + 3.86 * t9 ** (2. / 3)
                                          + .0027 * t9 ** -1 * exp(2.515E-3 * t9 ** -1))

    li_p = 1.096E9 * t9 ** (-2. / 3) * exp(-8.472 * t9 ** (-1. / 3)) - \
        4.830E8 * t9_star2 ** (5. / 6) * t9 ** (-3. / 2) * exp(-8.472 * t9_star2 ** (-1. / 3)) \
        + 1.06E10 * t9 ** (-3. / 2) * exp(-30.442 * t9 ** -1)
    """
    "These two reactions are not used in this assignment, but may be useful \
     at a later point, so I've kept them here."

    Berylliym_Boron = 	3.11E5*t9**(-2./3)*exp(-10.262*t9**(-1./3))+2.53E3\
                         *t9**(-3./2)*exp(-7.306*t9**-1)

    Nitrogen_Oxygen = 	4.9E7*t9**(-2./3)*exp(-15.228*t9**(-1./3)-\
                         .092*t9**2)*(1+.027*t9**(1./3)-.778*t9**(2./3)\
                         -.149*t9+.261*t9**(4./3)+.127*t9**(5./3)) + \
                         2.37E3*t9**(2./3)*exp(-3.011*t9**-1) + \
                         2.19E4*exp(-12.53*t9**-1)
    """
    rrs = [h_d, he3_p, he3_be, be_li, li_p]  # ,Berylliym_Boron,Nitrogen_Oxygen
    return rrs


# Computing reaction rates
"""
This segment computes the reaction rates for the different processes in the sun.
"""


def reaction_rates(t, rho, x, y3, y4, z_73li, z_74be):
    n = atom_numbers(rho, x, y3, y4, z_73li, z_74be)
    rrs = reactions(t)

    # Reaction rates
    """All of the following rates are calculated by the equation
    r_ik = n_i*n_k/(rho*(1+delta_ik)*lambda_ik)
    where i,k define the elements and delta_ik =1 if i=k else 0.
    """
    lambda_pp = rrs[0] / N_A * 1E-6  # [m^3/s]
    lambda_33 = rrs[1] / N_A * 1E-6  # [m^3/s]
    lambda_34 = rrs[2] / N_A * 1E-6  # [m^3/s]
    lambda_7e = rrs[3] / N_A * 1E-6  # [m^3/s]
    lambda_71prime = rrs[4] / N_A * 1E-6  # [m^3/s]

    r_pp = n[2] ** 2 / (rho * 2) * lambda_pp  # [kg-1*s-1]
    r_pd = r_pp  # [kg-1*s-1] Assume this reaction happens \
    # instantly such that it happens at the same \
    # rate as the elements it needs become available.
    # Thus it must be the same as the reaction that \
    # creates those elements (in this case Deuterium)
    r_33 = n[3] ** 2 / (rho * 2) * lambda_33  # [kg-1*s-1]
    r_34 = n[3] * n[4] / rho * lambda_34  # [kg-1*s-1]
    r_7e = n[5] * n[6] / rho * lambda_7e  # [kg-1*s-1]
    r_71prime = n[5] * n[2] / rho * lambda_71prime  # [kg-1*s-1]

    # This part makes sure that reaction rates which rely on other rates \
    # don't use elements that are not yet present.
    if 2. * r_33 + r_34 > r_pd:
        r_33 = 2. / 3 * r_pd
        r_34 = 1. / 3 * r_pd
    if r_7e > r_34:
        r_7e = r_34
    if r_71prime > r_7e:
        r_71prime = r_7e
    reaction_rate = r_pp, r_pd, r_33, r_34, r_7e, r_71prime
    return reaction_rate


# Solving the equation for the change in Luminosity
def energy(temp, rho, x, y3, y4, z_73li, z_74be):
    reaction_rate = reaction_rates(temp, rho, x, y3, y4, z_73li, z_74be)
    e_1 = reaction_rate[0] * (Q_pp + Q_dp)
    e_2 = reaction_rate[2] * Q_33
    e_3 = reaction_rate[3] * Q_34
    e_4 = reaction_rate[4] * Q_7e
    e_5 = reaction_rate[5] * Q_71prime
    e = e_1 + e_2 + e_3 + e_4 + e_5
    ppi_frac = e_2 + .69 * e_1
    ppii_frac = e_3 + e_4 + e_5 + .31 * e_1
    ppi_only = e_2
    ppii_only = e_3 + e_4 + e_5
    output = e, ppi_frac, ppii_frac, ppi_only, ppii_only
    return output


def opacity(t, rho, lines):
    # This function reads the opacity table
    # Convert r and T to the form used in opacity.txt
    r = (rho * 1E-3) / (t / 1E6)
    r_value = np.log10(r)
    t_value = np.log10(t)
    r_value = (round(2 * r_value) / 2.)
    # Print error message if using values outside table
    if r_value >= 1.5 or r_value <= -8.5:
        print 'Tried using r outside opacity table. Exiting'
        sys.exit(1)
    if t_value >= 8.8 or t_value <= 3.7:
        print 'Tried using T outside opacity table. Exiting'
        sys.exit(1)

    # Pick the right r value from the table
    for i in range(len(lines[0])):
        if r_value == lines[0][i]:
            r_i = i

    # Make list of T values found in table
    t_list = np.zeros(len(lines))
    for i in range(1, len(lines)):
        t_list[i] = lines[i][0]

    # Pick the right T value from the table
    t_i = min(range(len(t_list)), key=lambda j: abs(t_list[j] - t_value))

    # Pick the right kappa value with r and T
    log10kappa = lines[t_i][r_i]

    # Convert the value given to the form used in the equations
    kappa = .1 * 10. ** log10kappa
    return kappa


# Function to calculate xi
def xi_func(l_m, u, nabla_ad, nabla_rad):
    c1 = l_m ** 2 / u
    c2 = 1
    c3 = 4 * u / l_m ** 2
    c4 = nabla_ad - nabla_rad
    xi = np.roots([c1, c2, c3, c4])
    for element in xi:
        if np.imag(element) == 0:
            xi = np.real(element)
            break
    return xi
