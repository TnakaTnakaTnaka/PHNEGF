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
$ python phtherm.py --tran=(prefix).tran (--Tmin=0 --Tmax=1000 --dT=10)

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
BOLTHZ_CONST = 8.6173303e-5       # Boltzmann constant (ev / K)
HBAR = 6.582119514e-16  # Dirac constant (eV * s)
LIGHT_SPEED = 2.99792458e+10      # speed of light (cm / s)
EV_UNIT = 1.60217662e-19     # (J / eV)
DELTA = 1e-14           # parameter to avoid divergence in bose_function


def bose_function(x):
    return 1 / (np.exp(x + DELTA) - 1)


def dist_func(x):
    gx = np.exp(x) * (x * bose_function(x)) ** 2
    gx[0] = 1.0
    return gx


class PhThermCondCalculator:

    def __init__(self):
        self.__tran_file = ""
        self.__phtherm_file = ""
        self.__T_min = 0
        self.__T_max = 1000
        self.__T_width = 10


    def __set_parameters(self):
        options = parser.parse_args()
        if options.tran:
            self.__tran_file = options.tran
            prefix = self.__tran_file.split(".")[0]
            self.__phtherm_file = prefix + ".phtherm"
        else:
            print("tran file is not selected.")
            exit(1)
 
        if options.Tmin:
            self.__T_min = float(options.Tmin)
            print("The minimum tempreature : %3.1f K" % (self.__T_min))
            if self.__T_min < 0:
                print("Specified temperature is negative.")
                exit(1)
 
        if options.Tmax:
            self.__T_max = float(options.Tmax)
            print("The maximum tempreature : %3.1f K" % (self.__T_max))
            if self.__T_min > self.__T_max:
                print("Tmin is larger than Tmax. Check arguments.")
                exit(1)
 
        if options.dT:
            self.__T_width = float(options.dT)
            if self.__T_max == self.__T_min:
                self.__T_width = 0
            print("The tempreature width   : %3.1f K" % (self.__T_width))



    def __calc_phtherm(self, omega, tran, T):
        num_data = np.shape(omega)[0]
        beta = 1 / (BOLTHZ_CONST * T)
        omega_bar = 2.0 * math.pi * beta * HBAR * LIGHT_SPEED * omega  # dimensionless
        gx = dist_func(omega_bar) * tran  # Integrand
 
        delta_omega = omega[1] - omega[0]
 
        # Integration (Trapezoidal method)
        phtherm = 0.5 * (gx[0] + gx[num_data-1])
        for i in range(1, num_data-2):
            phtherm += gx[i]
 
        fac = EV_UNIT * BOLTHZ_CONST * LIGHT_SPEED * delta_omega
        phtherm *= fac  # unit : (W/K)
 
        return phtherm


    def start_calculation(self):

        self.__set_parameters() 

        # Load and initialize
        if self.__T_width != 0:
            step = int((self.__T_max - self.__T_min) / self.__T_width) + 1
        else:
            step = 1
        phtherm_data = np.zeros([step, 2]) # Initialize phonon thermal conductance data
        data = np.loadtxt(self.__tran_file)      # Load transmittance data file
        omega = data.T[0]                 # Frequency
        tran = data.T[1]                  # Transmittance
        
        # Loop for temperature
        for i in range(step):
            temperature = float(i * self.__T_width + self.__T_min)  # Temperature (K)
            if temperature > 0:
                # calculate thermal conductance
                phtherm = self.__calc_phtherm(omega, tran, temperature)
        
            else:
                phtherm = 0.0
        
            phtherm_data[i][0] = temperature
            phtherm_data[i][1] = phtherm
        
        # save kappa_phonon data
        np.savetxt(self.__phtherm_file, phtherm_data, delimiter='  ')
        print(self.__phtherm_file + " was generated.")



def main():

    # generate instance
    phtherm_calc = PhThermCondCalculator()
    phtherm_calc.start_calculation()
 
if __name__ == "__main__":
    main()
