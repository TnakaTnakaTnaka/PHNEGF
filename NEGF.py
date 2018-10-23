#!/usr/bin/env python
#
# NEGF.py is a sclipt to obtain the phonon transmission function
# for each k-point by using the hessian file in ALAMODE.
#
#
# Copyright (c) 2018 Yuto Tanaka
#

"""
--- How to use ---
$ python NEGF.py --negf=(prefix)_negf.in --hessian=(prefix).hessian

"""

import argparse
import time
import numpy as np
import numpy.linalg as LA
import mod_dymat as dymat

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage=usage)
parser.add_argument("--negf", help="negf file")
parser.add_argument("--hessian", help="hessian file")


def surface_green(nat, OMG_s, D_s):

    N_atom = 3 * nat

    ep_s = OMG_s - (D_s[:N_atom, :N_atom])
    ep = OMG_s - (D_s[:N_atom, :N_atom])
    alpha = (D_s[:N_atom, N_atom:])
    beta = (D_s[N_atom:, :N_atom])
    G_0 = LA.inv(ep_s)
    norm_G = 1.0
    while norm_G > criterion:
        ep_inv = LA.inv(ep)
        a = np.dot(ep_inv, alpha)
        b = np.dot(ep_inv, beta)
        ep_s -= np.dot(alpha, b)
        ep -= np.dot(beta, a) + np.dot(alpha, b)
        alpha = np.dot(alpha, a)
        beta = np.dot(beta, b)

        G = LA.inv(ep_s)
        norm_G = LA.norm(G-G_0)
        G_0 = G

    return G


def transmission(i, nat, D_c, D_s, D_cl, D_cr):

    omega = i * grid
    omega2 = (omega)**2 + 1e-10
    OMG = omega2 * (1 + delta*1j)

    # Frequency
    OMG_c = OMG * np.identity(9 * nat, dtype=np.complex128)
    OMG_s = OMG * np.identity(3 * nat, dtype=np.complex128)

    # coupling term in dynamical matrix
    D_lc = np.conjugate(D_cl.T)
    D_rc = np.conjugate(D_cr.T)

    # Surface Green's function
    G_l = surface_green(nat, OMG_s, D_s[::-1, ::-1])[::-1, ::-1]
    G_r = surface_green(nat, OMG_s, D_s)

    # Self energy
    Self_l = np.dot(D_cl, np.dot(G_l, D_lc))
    Self_r = np.dot(D_cr, np.dot(G_r, D_rc))

    # Green's function in the scattering ragion
    G_c = LA.inv(OMG_c - D_c - (Self_l + Self_r))
    G_c_her = np.conjugate(G_c.T)

    # Gamma
    Gamma_l = (Self_l - np.conjugate(Self_l.T)) * 1j
    Gamma_r = (Self_r - np.conjugate(Self_r.T)) * 1j

    return omega * cm, np.trace(np.dot(Gamma_l, np.dot(G_c, np.dot(Gamma_r, G_c_her)))).real


def generate_qmesh(kpoint, tran_direct):

    num_q = kpoint[0] * kpoint[1] * kpoint[2]
    q = np.zeros([num_q, 3])

    var_idx = [i for i, x in enumerate(tran_direct) if x == 0]
    fix_idx = [i for i, x in enumerate(tran_direct) if x == 1][0]
    # fix_idx = tran_direct.index(1)
    bz = [[], [], []]
    bz[fix_idx].append(0.0)

    for i in var_idx:
        dq = 1 / (kpoint[i] + 1)
        qx = -0.5 + dq

        while qx < 0.5 - 1e-3:
            bz[i].append(qx)
            qx += dq

    q_count = 0
    for i in range(kpoint[0]):
        for j in range(kpoint[1]):
            for k in range(kpoint[2]):
                q[q_count] = np.array([bz[0][i], bz[1][j], bz[2][k]])
                q_count += 1

    return q, var_idx


def get_qpoint(qmesh, revec):
    return qmesh[0] * revec[0] + qmesh[1] * revec[1] + qmesh[2] * revec[2]


def main():

    # grobalization
    global cm
    global delta
    global criterion
    global grid

    start = time.time()

    options = parser.parse_args()
    if options.negf:
        negf_file = options.negf
    else:
        print("negf file is not selected.")

    if options.hessian:
        hessian_file = options.hessian
    else:
        print("hessian file is not selected.")

    prefix = negf_file.split('.')[0]

    # read negf file
    x_bohr, k_atom, nat, mass, lavec, univec, revec, tran_direct, \
            kpoint, cutoff, delta, freq_max, criterion, step \
            = dymat.read_negf(negf_file)

    # atoms in unitcell  atom_uc = [1, ..., nat_unitcell]
    atom_uc = dymat.atom_in_unitcell(x_bohr, univec, nat)
    nat_uc = len(atom_uc)  # number of atoms in unit cell
    
    # supercell infomation
    lmn = dymat.supercell(lavec, univec)

    # shift parameter
    dymat.make_shift_list(lmn)

    # considerable atom pairs for fcs
    pairs = dymat.generate_pairs(atom_uc, x_bohr, lavec, univec, nat, cutoff)

    # mapping equivalant atom in unit cell
    map_uc = dymat.mapping(x_bohr, univec, atom_uc, nat, lmn)

    # atomic mass in uni tcell
    mass_uc = dymat.mass_in_unitcell(mass, k_atom, atom_uc)

    # store fcs matrix
    fcs = dymat.store_all_fcs(hessian_file, atom_uc, nat_uc, pairs, map_uc, mass_uc, tran_direct)

    # obtain k-point in 1st BZ and transport direction index
    qmesh, var_idx = generate_qmesh(kpoint, tran_direct)

    q_count = 0
    cm = 3634.87331806918 # convert to cm^{-1}
    freq_max /= cm
    grid = float(freq_max) / step

    for i in range(kpoint[var_idx[0]]):
        for j in range(kpoint[var_idx[1]]):
            outfile = prefix + ".tran" + str(i) + "_" + str(j)
            print(outfile)

            tran_data = np.zeros([step, 2])

            q = get_qpoint(qmesh[q_count], revec)
            q_count += 1

            #Dynamical matrix
            D_c, D_s, D_cl, D_cr = dymat.generate_dynamical_matrix( \
                    fcs, q, nat_uc, univec, tran_direct)

            #for loop frequency
            for s in range(step):
                omega, Tran = transmission(s, nat_uc, D_c, D_s, D_cl, D_cr)
                tran_data[s][0] = omega
                tran_data[s][1] = Tran

                np.savetxt(outfile, tran_data, delimiter='  ')

    print(time.time() - start, "seconds")

if __name__ == "__main__":
    main()
