#!/usr/bin/env python
#
# kappa.py
#
# script to calculate thermal conductance from the result of transmittance
#
#
# Copyright (c) 2018 Yuto Tanaka
#

"""
--- How to use ---
$ python kappa.py --tran=(prefix).tran (--Tmin=0 --Tmax=1000 --dT=10)

--- default parameter---
Tmin = 0 K
Tmax = 1000 K
dT   = 10 K
"""

import argparse
import math
import numpy as np

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage=usage)
parser.add_argument("tran", help="tran file")
parser.add_argument("--Tmin", action="store", default="0",
                    help="print the minimum temperature you want to calculate.")
parser.add_argument("--Tmax", action="store", default="1000",
                    help="print the maximum temperature you want to calculate.")
parser.add_argument("--dT", action="store", default="10",
                    help="print the width of temperature. Default is 10 K")


# parameters
kb = 8.6173303e-5       # Boltzmann constant (ev / K)
hbar = 6.582119514e-16  # Dirac constant (eV * s)
c = 2.99792458e+10      # speed of light (cm / s)
eV = 1.60217662e-19     # (J / eV)
delta = 1e-14           # parameter to avoid divergence in bose_function


def bose_function(x):
    return 1 / (np.exp(x + delta) - 1)


def dist_func(x):
    gx = np.exp(x) * (x * bose_function(x)) ** 2
    gx[0] = 1.0
    return gx


def calc_kappa(omega, tran, T):
    num_data = np.shape(omega)[0]
    beta = 1 / (kb * T)
    omega_bar = 2.0 * math.pi * beta * hbar * c * omega  # dimensionless
    gx = dist_func(omega_bar) * tran  # Integrand

    domega = omega[1] - omega[0]

    # Integration (Trapezoidal method)
    k_p = 0.5 * (gx[0] + gx[num_data-1])
    for i in range(1, num_data-2):
        k_p += gx[i]

    fac = eV * kb * c * domega
    k_p *= fac  # unit : (W/K)

    return k_p


def main():

    # Set args
    options = parser.parse_args()
    if options.tran:
        tran_file = options.tran
        prefix = tran_file.split(".")[0]
        kappa_file = prefix + ".kl"
    else:
        print("tran file is not selected.")
        exit(1)

    if options.Tmin:
        T_min = float(options.Tmin)
        print("The minimum tempreature : %3.1f K" % (T_min))

    if options.Tmax:
        T_max = float(options.Tmax)
        print("The maximum tempreature : %3.1f K" % (T_max))

    if options.dT:
        T_width = float(options.dT)
        if T_max == T_min:
            T_width = 0
        print("The tempreature width   : %3.1f K" % (T_width))

    # Load and initialize
    if T_width != 0:
        step = int((T_max - T_min) / T_width) + 1
    else:
        step = 1
    kappa_data = np.zeros([step, 2])  # Initialize thermal conductance data
    data = np.loadtxt(tran_file)      # Load transmittance data file
    omega = data.T[0]                 # Frequency
    tran = data.T[1]                  # Transmittance

    # Loop for temperature
    for i in range(step):
        T = float(i * T_width + T_min)  # Temperature (K)
        if T > 0:
            # calculate thermal conductance
            kappa_p = calc_kappa(omega, tran, T)

        else:
            kappa_p = 0.0

        kappa_data[i][0] = T
        kappa_data[i][1] = kappa_p

    # save kappa_phonon data
    np.savetxt(kappa_file, kappa_data, delimiter='  ')
    print(kappa_file + " was generated.")


if __name__ == "__main__":
    main()
