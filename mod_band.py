#!/usr/bin/env python
#
# mod_band.py
#
# Copyright (c) 2018 Yuto Tanaka
#

import cmath
import numpy as np
import numpy.linalg as LA


cm = 3634.87331806918 # convert to cm^{-1}

def calc_band(nat_uc, univec, revec, D_c):
    N = 3 * nat_uc
    g = 0.0
    X = 0.5
    step = 200
    dq = float( (X-g) / step )
    qx = 0.0

    dymat = np.zeros([3, N, N], dtype=np.complex128)
    dymat[-1] = D_c[N:2*N, :N]
    dymat[0]  = D_c[N:2*N, N:2*N]
    dymat[1]  = D_c[N:2*N, 2*N:]
    vec_l = -univec[0] 
    vec_r =  univec[0]

    file_band = 'prefix.bands'
    f_out = open(file_band, 'w')
    f_out.write('#  G X\n')
    f_out.write('#  0.00000 %6.5f\n' %(X * revec[0][0]))
    f_out.write('#  k-axis, Eigenvalues [cm^-1]\n')
    q_p = np.array([0, 0.0, 0.0])
    while qx <= X + dq:
        D = np.zeros([N, N], dtype=np.complex128)
        q = qx * revec[0] + q_p
 
        theta = np.zeros([3])
        theta[-1] = np.dot(q, vec_l)
        theta[1]  = np.dot(q, vec_r)

        for i in [-1, 0, 1]:
            D += dymat[i] * cmath.exp(theta[i] * 1j)

        eig_val = LA.eigvalsh(D)
        eig_val[eig_val < 0] = 0
        # eig_val[abs(eig_val) < 1e-7] = 0.0
        if np.all(eig_val >= 0.0):
            eig = np.sqrt(eig_val).real * cm
 
        else:
            print("negative frequency exists.")
            exit(1)

        f_out.write('%f  ' % (q[0]))
        f_out.write(' '.join(map(str, eig)) + '\n')

        qx += dq

    f_out.close()

