#!/usr/bin/env python
#
# mod_dymat.py
#
# module for generating dynamical matrix to calculate phonon
# transmission function by using NEGF method.
#
# Copyright (c) 2018 Yuto Tanaka
#

import math
import cmath
import numpy as np
import numpy.linalg as LA


def reciprocal(lavec):

    revec = np.zeros([3, 3])
    vol = np.dot(lavec[0], (np.cross(lavec[1], lavec[2])))

    revec[0] = np.cross(lavec[1], lavec[2])
    revec[1] = np.cross(lavec[2], lavec[0])
    revec[2] = np.cross(lavec[0], lavec[1])

    revec *= 2 * math.pi / vol

    return revec


def store_vec(factor, vec, vec_count, ss):

    if factor == 0:
        factor = float(ss[0])

    else:
        for i in range(3):
            vec[vec_count][i] = factor * float(ss[i])

        vec_count += 1

    return factor, vec, vec_count


def validation(tran_direct, kpoint):

    val = True

    var_idx = [i for i, x in enumerate(tran_direct) if x == 0]
    fix_idx = [i for i, x in enumerate(tran_direct) if x == 1]
    if len(var_idx) != 2 or len(fix_idx) != 1:
        print("Transport direction is not selected properly.")
        print("Please check direction field in negf input file.")
        val = False

    kp_wrong = [i for i, x in enumerate(kpoint) if x < 0]
    if len(kp_wrong) > 0:
        print("k point must be positive integer.")
        print("Please check kpoint field in negf input file.")
        val = False

    if not val:
        exit(1)


def read_negf(negf_file):

    target = ['nat', 'nkd', 'mass', '&cell', '&unit_cell',
              '&direction', '&kpoint', '&position', '&cutoff']
    negf_target = ['imag_delta', 'freq_max', 'criterion', 'freq_div']
    lavec_frag = 0
    univec_frag = 0
    pos_frag = 0
    direct_frag = 0
    kpoint_frag = 0
    cutoff_frag = 0

    lavec = np.zeros([3, 3])  # supercell lattice vector (initialize)
    univec = np.zeros([3, 3])  # unit vector (initialize)

    # default NEGF parameters
    imag_delta = 1e-6
    criterion = 1e-6
    freq_max = 1000
    step = 100

    f_negf = open(negf_file, 'r')  # open output file
    for line in f_negf:
        ss = line.strip().split()
        if len(ss) == 0:  # empty line
            continue

        else:
            # supercell lattice vector
            if lavec_frag == 1:
                factor, lavec, vec_count = store_vec(
                    factor, lavec, vec_count, ss)

                if vec_count == 3:
                    lavec_frag = 0
                    conv = LA.inv(LA.inv(lavec).T)  # unit converter

            # unit cell vector and reciprocal vector
            elif univec_frag == 1:
                factor, univec, vec_count = store_vec(
                    factor, univec, vec_count, ss)

                if vec_count == 3:
                    univec_frag = 0
                    revec = reciprocal(univec)  # reciprocal vector

            # transport direction
            elif direct_frag == 1:
                tran_direct = [int(ss[0]), int(ss[1]), int(ss[2])]
                direct_frag = 0

            # kpoint
            elif kpoint_frag == 1:
                if kp_mode < 0:
                    kp_mode = int(ss[0])
                    if kp_mode != 3:
                        print(
                            'If you calculate dynamical matrix, KPMODE should be 3.')
                        exit(1)

                else:
                    kpoint = [int(ss[0]), int(ss[1]), int(ss[2])]
                    kpoint_frag = 0

            # coordinates (bohr unit)
            elif pos_frag == 1:
                k_atom[atom_count] = int(ss[0])
                for i in range(3):
                    x_bohr[atom_count][i] = float(ss[i+1])

                atom_count += 1
                if atom_count == nat:
                    pos_frag = 0

            elif cutoff_frag == 1:
                cutoff = ss[1]
                cutoff_frag = 0

            # search target and build frag
            if ss[0].lower() == target[0]:
                nat = int(ss[2])  # number of atom
                atom_count = 0
                # kind of atom (initialize)
                k_atom = np.zeros([nat], dtype=np.int64)
                x_bohr = np.zeros([nat, 3])  # atomic coordinate (initialize)

            elif ss[0].lower() == target[1]:
                nkd = int(ss[2])  # number of kind of atom
                mass = np.zeros([nkd])  # atomic mass (initialize)

            elif ss[0].lower() == target[2]:
                for i in range(nkd):
                    mass[i] = ss[i+2]  # atomic mass

            elif ss[0].lower() == target[3]:
                factor = 0
                vec_count = 0
                lavec_frag = 1

            elif ss[0].lower() == target[4]:
                factor = 0
                vec_count = 0
                univec_frag = 1

            elif ss[0].lower() == target[5]:
                direct_frag = 1
                tran_direct = []

            elif ss[0].lower() == target[6]:
                kpoint_frag = 1
                kpoint = []
                kp_mode = -1

            elif ss[0].lower() == target[7]:
                pos_frag = 1

            elif ss[0].lower() == target[8]:
                cutoff_frag = 1

            # NEGF options
            elif ss[0].lower() == negf_target[0]:
                imag_delta = float(ss[2])

            elif ss[0].lower() == negf_target[1]:
                freq_max = float(ss[2])

            elif ss[0].lower() == negf_target[2]:
                criterion = float(ss[2])

            elif ss[0].lower() == negf_target[3]:
                step = int(ss[2])

    f_negf.close()

    # convert coordinate unit frac to bohr
    for i in range(nat):
        x_bohr[i] = np.dot(conv, x_bohr[i])

    # set cutoff radius
    if cutoff.lower() == 'none':
        cutoff = univec[0][0]
    else:
        cutoff = float(cutoff)

    validation(tran_direct, kpoint)

    return x_bohr, k_atom, nat, mass, lavec, univec, revec, \
        tran_direct, kpoint, cutoff, imag_delta, freq_max, \
        criterion, step


def supercell(lavec, univec):
    lmn = [int(lavec[i][i] / univec[i][i]) for i in range(3)]

    return lmn


def make_shift_list(lmn):

    shift_max = [1, 1, 1]  # shift_max = [max_x, max_y, max_z]
    for i in range(3):
        if lmn[i] > 3:
            shift_max[i] = int(lmn[i] * 0.5)

    global xrng_u, yrng_u, zrng_u
    global lavec_shift, univec_shift

    def make_rng(m):
        rng_l = [i for i in range(-m, 1)]
        rng_u = [i for i in range(-m, m+1)]
        return rng_l, rng_u

    xrng_l, xrng_u = make_rng(shift_max[0])
    yrng_l, yrng_u = make_rng(shift_max[1])
    zrng_l, zrng_u = make_rng(shift_max[2])

    lavec_shift = [[i, j, k] for i in xrng_l for j in yrng_l for k in zrng_l]
    univec_shift = [[i, j, k] for i in xrng_u for j in yrng_u for k in zrng_u]


def in_unitcell(x, vec):
    d = 1e-7
    return 0 <= x[0] < vec[0]-d and 0 <= x[1] < vec[1]-d and 0 <= x[2] < vec[2]-d


def calc_shift(shift, vec):
    return shift[0] * vec[0] + shift[1] * vec[1] + shift[2] * vec[2]


def atom_in_unitcell(x_bohr, univec, nat):
    unit_len = np.sum(univec, axis=0)
    atom_uc = [i + 1 for i in range(nat) if in_unitcell(x_bohr[i], unit_len)]

    return atom_uc


def check_symmetry(lavec, univec, x_bohr, p, q, cutoff):
    shift = np.zeros([3])
    distance = []
    symmetry = []
    unit_len = np.sum(univec, axis=0)

    Nshift = len(lavec_shift)
    for i in range(Nshift):
        shift = calc_shift(lavec_shift[i], lavec)
        dist = LA.norm(x_bohr[p-1] - (x_bohr[q-1] + shift))
        distance.append(dist)

    if min(distance) > cutoff:
        pass

    else:
        index = [i for i, x in enumerate(
            distance) if abs(x - min(distance)) < 1e-5]
        Nshift = len(univec_shift)
        for i in index:
            shift = calc_shift(lavec_shift[i], lavec)
            x_la = x_bohr[q-1] + shift
            for j in range(Nshift):
                shift = calc_shift(univec_shift[j], univec)
                x_uni = x_la - shift

                if in_unitcell(x_uni, unit_len):
                    symmetry.append(j)
                    break

    return symmetry


def generate_pairs(atom_uc, x_bohr, lavec, univec, nat, cutoff):

    pairs = {}

    for p in atom_uc:
        pairs[p] = []
        for q in range(p, nat + 1):
            sym = check_symmetry(lavec, univec, x_bohr, p, q, cutoff)
            pairs[p].append([q, sym])

    f_log = open("pairs.log", "w")
    f_log.write("atom number in the unit cell\n")
    f_log.write(str(atom_uc) + "\n\n")

    f_log.write("univec_shift R = l*a + m*b + n*c\n")
    num_shift = len(univec_shift)
    for i in range(num_shift):
        f_log.write(str(i+1) + " " + str(univec_shift[i]) + "\n")
    f_log.write("\n")

    for p in atom_uc:
        f_log.write("symmetry pair of atom %d\n" % (p))
        for q in range(len(pairs[p])):
            f_log.write(str(pairs[p][q]) + "\n")
        f_log.write("\n")
    f_log.close()

    return pairs


def mass_in_unitcell(mass, k_atom, atom_uc):
    mass_uc = [mass[k_atom[i-1] - 1] for i in atom_uc]

    return mass_uc


def mapping(x_bohr, univec, atom_uc, nat, lmn):
    map_uc = []
    unit_len = np.sum(univec, axis=0)
    for x in range(nat):
        break_frag = 0
        for i in range(lmn[0]):
            for j in range(lmn[1]):
                for k in range(lmn[2]):
                    shift = calc_shift([i, j, k], univec)
                    x_uni = x_bohr[x] - shift
                    if in_unitcell(x_uni, unit_len):
                        break_frag = 1
                        break  # break for loop k

                # break for loop j
                if break_frag == 1:
                    break

            # break for loop j
            if break_frag == 1:
                break

        for p in atom_uc:
            r = LA.norm(x_uni - x_bohr[p-1])
            if r < 1e-3:
                map_uc.append(p)
                break

        if len(map_uc) != x + 1:
            print("mapping atom error.")
            exit(1)

    return map_uc


def store_all_fcs(hessian_file, atom_uc, nat_uc, pairs, map_uc, mass_uc):

    fcs = np.zeros([len(xrng_u), len(yrng_u), len(zrng_u),
                    3*nat_uc, 3*nat_uc], dtype=np.complex128)

    num_uniq_pair = {i: len(pairs[i]) for i in pairs.keys()}

    f_in = open(hessian_file, 'r')
    label = f_in.readline()   # Do not comment out
    for line in f_in:
        ss = line.strip().split()

        atom1 = int(ss[0])
        xyz1 = int(ss[1])
        atom2 = int(ss[2])
        xyz2 = int(ss[3])
        fc2 = float(ss[4])

        if fc2 == 0.0 or atom1 not in atom_uc or atom1 > atom2:
            continue

        else:
            idx1_uc = atom_uc.index(map_uc[atom1 - 1])
            idx2_uc = atom_uc.index(map_uc[atom2 - 1])
            idx1 = idx1_uc + (xyz1 - 1) * nat_uc
            idx2 = idx2_uc + (xyz2 - 1) * nat_uc

            for i in range(num_uniq_pair[atom1]):
                if pairs[atom1][i][0] == atom2:
                    p_idx = i
                    break

            if pairs[atom1][p_idx][1] != []:
                dev = len(pairs[atom1][p_idx][1]) * \
                    math.sqrt(mass_uc[idx1_uc] * mass_uc[idx2_uc])
                fc2 /= dev

            else:
                continue

            for i in pairs[atom1][p_idx][1]:
                idx_x = univec_shift[i][0]
                idx_y = univec_shift[i][1]
                idx_z = univec_shift[i][2]
                fcs[+idx_x][+idx_y][+idx_z][idx1][idx2] = fc2
                fcs[-idx_x][-idx_y][-idx_z][idx2][idx1] = fc2

    f_in.close()

    return fcs


def generate_dynamical_matrix(fcs, q, nat, univec, tran_direct):
    N = 3 * nat
    dymat = np.zeros([3, N, N], dtype=np.complex128)

    # initial values
    D_c = np.zeros([3*N, 3*N], dtype=np.complex128)
    D_s = np.zeros([2*N, 2*N], dtype=np.complex128)
    D_cl = np.zeros([3*N, N], dtype=np.complex128)
    D_cr = np.zeros([3*N, N], dtype=np.complex128)

    # transport along x direction
    if tran_direct[0] == 1:
        for k in xrng_u:
            for l in yrng_u:
                for m in zrng_u:
                    r_vec = l * univec[1] + m * univec[2]
                    theta = np.dot(q, r_vec)
                    dymat[k] += fcs[k][l][m] * cmath.exp(theta * 1j)

    # transport along y direction
    elif tran_direct[1] == 1:
        for l in yrng_u:
            for m in zrng_u:
                for k in xrng_u:
                    r_vec = m * univec[2] + k * univec[0]
                    theta = np.dot(q, r_vec)
                    dymat[l] += fcs[k][l][m] * cmath.exp(theta * 1j)

    # transport along z direction
    elif tran_direct[2] == 1:
        for m in zrng_u:
            for k in xrng_u:
                for l in yrng_u:
                    r_vec = k * univec[0] + l * univec[1]
                    theta = np.dot(q, r_vec)
                    dymat[m] += fcs[k][l][m] * cmath.exp(theta * 1j)

    else:
        print('Transport direction is not selected.')
        exit(1)

    for i in range(3):
        index = i * N
        D_c[index:index + N, index:index + N] = dymat[0]

    D_c[:N, N:2*N] = np.conjugate(dymat[-1].T)
    D_c[N:2*N, :N] = dymat[-1]
    D_c[N:2*N, 2*N:] = dymat[+1]
    D_c[2*N:, N:2*N] = np.conjugate(dymat[+1].T)

    D_s = D_c[:2*N, :2*N]
    D_cl[:N, :] = dymat[-1]
    D_cr[2*N:, :] = dymat[1]

    return D_c, D_s, D_cl, D_cr
