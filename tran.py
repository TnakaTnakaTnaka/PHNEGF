#!/usr/bin/env python
#
# tran.py
#
# script to calculate transmittance averaged in 1st Brillouin zone
# from the data obtained by NEGF.py.
#
#
# Copyright (c) 2018 Yuto Tanaka
#

import os
import sys
import fnmatch
import numpy as np


class TransmissionCalculator:

    def __init__(self):
        pass

    def __data_column(self):
        # find k-resolved transmission file
        flag = False;
        for file_name in os.listdir(os.getcwd()):
            if fnmatch.fnmatch(file_name, '*.tran0_0'):
                prefix = file_name.split('.')[0]
                data = np.loadtxt(file_name)
                frequency = data.T[0]
                num_data = int(np.shape(data)[0])
                flag = True
 
        # If the *.tran0_0 file is not exist, exit this program
        if False == flag:
            print("*.tran0_0 file is not exist.")
            sys.exit()
 
        return num_data, frequency, prefix


    def __integration(self, num_data):
        # integrate transmission function in 1st BZ
        tran = np.zeros(num_data)
        kpoint = 0
        for file_name in os.listdir('.'):
            if fnmatch.fnmatch(file_name, '*.tran*_*'):
                kpoint += 1  # kpoint countor
                data = np.loadtxt(file_name)
                tran += data.T[1]
 
        tran /= kpoint
        print(kpoint, "file(s) is loaded.")
 
        return tran


    def start_calculation(self):
        # Set transmission array
        num_data, frequency, prefix = self.__data_column()
        tran = np.zeros([num_data, 2])

        tran.T[0] = frequency
        tran.T[1] = self.__integration(num_data)

        # Save transmission file
        outfile = prefix + '.tran'
        np.savetxt(outfile, tran, delimiter='  ')
        print("%s was generated." % (outfile))



def main():

    # Create Instance
    tran_calc = TransmissionCalculator()

    # Start transmission calculation
    tran_calc.start_calculation();

if __name__ == "__main__":
    main()
