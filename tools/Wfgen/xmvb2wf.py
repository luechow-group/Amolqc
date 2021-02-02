#!/usr/bin/env python

"""
 Copyright (C) 2017-2018 Leonard Reuter

SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
from .Utils import *


def xmvb_in(input_name, info_name, basis, wf_type):
    xmo_filename = input_name + '.xmo'
    info_filename = info_name
    orbital_format = 'gau'

    bohr = 0.52917721067  # 2014 CODATA

    # opening files
    try:
        info_file = open(info_filename)
    except FileNotFoundError:
        not_found(info_filename)

    try:
        xmo_file = open(xmo_filename)
    except FileNotFoundError:
        not_found(xmo_filename)

    # reading geometry and number of basis functions from INFO file
    line = info_file.readline()
    words = line.split()
    number_atoms = int(words[0])
    number_basisfunctions = int(words[1])

    charge = 0
    number_inactive_electrons = 0

    atoms, charge, line = read_atoms(bohr, charge, info_file, line, number_atoms)

    info_file.close()

    # reading localization from .xmo file
    bovb = False
    while '$ctrl' not in line:
        line = xmo_file.readline()
        if not line:
            sys.exit('Error: no ctrl section found in .xmo file.')
    line = xmo_file.readline()
    while '$end' not in line:
        if 'bovb' in line.lower():
            bovb = True
        line = xmo_file.readline()

    while '$frag' not in line:
        line = xmo_file.readline()
        if not line:
            sys.exit('Error: no frag section found in .xmo file.')
    for i in range(2):
        line = xmo_file.readline()
    fragments = []
    while '$end' not in line:
        fragments.append(line.split('\n')[0])
        line = xmo_file.readline()

    total_symmetry = 'VB'

    while '$orb' not in line:
        line = xmo_file.readline()
        if not line:
            sys.exit('Error: no orb section found in .xmo file.')
    for i in range(2):
        line = xmo_file.readline()

    locality_index = [0] * len(fragments)
    orbitals = []
    while '$end' not in line:
        if line.split() and line.split()[0][0] != '#':
            fragment = int(line.split()[0])
            locality_index[fragment - 1] += 1
            orbital = Orbital()
            orbital.symmetry = str(locality_index[fragment - 1]) + '.' + str(fragment)
            orbitals.append(orbital)
        line = xmo_file.readline()

    if bovb:
        while 'Breathing Orbitals' not in line:
            line = xmo_file.readline()
            if not line:
                sys.exit('Error: BOVB information not found in .xmo file.')
        line = xmo_file.readline().split(':')[-1]
        while 'Reading' not in line:
            orbital = Orbital()
            fragment = int(orbitals[int(line.split('->')[0]) - 1].symmetry.split('.')[1])
            locality_index[fragment - 1] += 1
            orbital.symmetry = str(locality_index[fragment - 1]) + '.' + str(fragment)
            orbitals.append(orbital)
            line = xmo_file.readline().split(':')[-1]

    # reading charge, number and type of electrons from .xmo file
    while 'The following structures are used in calculation:' not in line:
        line = xmo_file.readline()
        if not line:
            sys.exit('Error: no used structures found in .xmo file.')
    for i in range (2):
        line = xmo_file.readline()
    words = line.split()
    if ':' in words[2]:
        # if the actual number of inactives is 2, this variable is 0
        number_inactive_electrons = 2 * int(words[2].split(':')[1])
        number_electrons = number_inactive_electrons + len(words) - 3
    else:
        number_electrons = len(words)-2
    line = xmo_file.readline()
    words = line.split()
    while not (not words or words[1] == '*****'):
        number_electrons += len(words)
        line = xmo_file.readline()
        words = line.split()
    charge -= number_electrons

    # read multiplicity and prepare reading determinants
    active_alpha_electrons = 0
    active_beta_electrons = 0
    determinant_lines = 0
    multiplicity = 0
    csfs = []
    # in the xmvb 3.0 output, the orbital format is slightly different
    xmvb3 = 1
    if wf_type != 'orb':
        while 'COEFFICIENTS OF DETERMINANTS' not in line:
            line = xmo_file.readline()
            if not line:
                sys.exit('Error: determinants not found in.xmo file')
        if 'UNNORMALIZED' not in line:
            xmvb3 = 0
            print("Warning: in XMVB 2 or earlier, unnormalized coefficients of determinants are not available."
                  "the read in wave function is thus not equal to the XMVB wave function.")
        for i in range (2):
            line = xmo_file.readline()
        while '******' not in line:
            words = line.split()
            for i in range(len(words)):
                if words[i] == 'a':
                    active_alpha_electrons += 1
                elif words[i] == 'b':
                    active_beta_electrons += 1
            line = xmo_file.readline()
            determinant_lines += 1
        multiplicity = active_alpha_electrons - active_beta_electrons + 1
        if multiplicity != 1 and wf_type == 'csf':
            print("The multiplicity is '"+str(multiplicity)+"'.")
            print('Warning: Do not trust csf generation for non-singlet systems.')
        number_alpha_electrons = (number_electrons + multiplicity - 1) // 2
        number_beta_electrons = number_electrons - number_alpha_electrons

        csfs = []

        # reading determinants
        csf = Csf()
        alt_csf = Csf()
        alpha_lines = determinant_lines // 2
        beta_lines = determinant_lines - alpha_lines - 1
        words = line.split()
        while words:
            determinant = Determinant()
            # reading determinant coefficient
            determinant.coefficient = float(words[1])
            for i in range(3):
                del words[0]
            for i, sign in zip([alpha_lines, beta_lines], [1, -1]):
                if number_inactive_electrons != 0:
                    for j in range(number_inactive_electrons // 2):
                        determinant.orbital_list.append(sign * (j + 1))
                    del words[0]
                for j in range(i):
                    for k in range(len(words)):
                        determinant.orbital_list.append(sign * int(words[k]))
                    line = xmo_file.readline()
                    words = line.split()

            alt_determinant = Determinant()
            alt_determinant.orbital_list = list(determinant.orbital_list)
            alt_determinant.coefficient = float(determinant.coefficient)
            alt_determinant.sort(change_sign=False)
            alt_csf.determinants.append(alt_determinant)

            determinant.sort()
            csf.determinants.append(determinant)
        csfs.append(csf)

    if wf_type == 'csf':
        det_csf = alt_csf
        csfs = []
        # reading rumer structures
        xmo_file.seek(0)
        while 'COEFFICIENTS OF STRUCTURES' not in line:
            line = xmo_file.readline()
            if not line:
                sys.exit('Error: coefficients of structures not found in .xmo file')
        for i in range(2):
            line = xmo_file.readline()
        number_csfs = 0
        structures = []
        electron_list = []
        words = line.split()

        while words:
            determinant = Determinant()
            if abs(float(words[1])) > 0:
                determinant.coefficient = float(words[1])
            else:
                determinant.coefficient = 0.0000001
            for i in range(3):
                del words[0]
            if number_inactive_electrons != 0:
                for i in range(number_inactive_electrons // 2):
                    determinant.orbital_list.append(i + 1)
                del words[0]
            while words and (len(words) < 3 or words[2] != '******'):
                for i in range(len(words)):
                    electron_list.append(words[i])
                line = xmo_file.readline()
                words = line.split()
            for i in range(active_alpha_electrons):
                if i < active_beta_electrons:
                    j = i
                determinant.orbital_list.append(int(electron_list[j]))
                del electron_list[j]
            for i in range(number_inactive_electrons // 2):
                determinant.orbital_list.append(-(i + 1))
            for i in range(len(electron_list)):
                determinant.orbital_list.append(-1 * int(electron_list[i]))
            electron_list = []
            number_csfs += 1
            structures.append(determinant)
        # building csfs from rumer structures
        for i in range(number_csfs):
            number_covalent_bonds = 0
            # building determinants of csfs
            csf = Csf()

            singly_occupied = []
            csf.coefficient = structures[i].coefficient
            for j in range(active_beta_electrons):
                if structures[i].orbital_list[number_inactive_electrons // 2 + j] \
                        != -1 * structures[i].orbital_list[number_inactive_electrons // 2 + j + number_alpha_electrons]:
                    number_covalent_bonds += 1
                    singly_occupied.append(number_inactive_electrons // 2 + j)

            csf.number_determinants = 2 ** number_covalent_bonds

            for j in range(2 ** number_covalent_bonds):
                determinant = Determinant()
                determinant.orbital_list = list(structures[i].orbital_list)
                determinant.coefficient = 2 ** (-0.5 * number_covalent_bonds)
                permutation_list = number2binarylist(j)
                for k in range(len(permutation_list)):
                    if permutation_list[k] == 1:
                        alpha_position = singly_occupied[k]
                        beta_position = singly_occupied[k] + number_alpha_electrons
                        determinant.orbital_list[alpha_position], determinant.orbital_list[beta_position] \
                            = -1 * determinant.orbital_list[beta_position], -1 * determinant.orbital_list[
                            alpha_position]
                determinant.sort()
                csf.determinants.append(determinant)

            # occupation is atm not used
            csf.occupation = [0] * len(orbitals)
            for j in range(number_beta_electrons):
                csf.occupation[structures[i].orbital_list[j] - 1] =\
                    -1 * structures[i].orbital_list[j + number_alpha_electrons]
                csf.occupation[-1 * structures[i].orbital_list[j + number_alpha_electrons] - 1] =\
                    structures[i].orbital_list[j]
            for j in range(multiplicity - 1):
                csf.occupation[structures[i].orbital_list[j + number_beta_electrons] - 1] = -1

            csfs.append(csf)

        # determining determinant coefficients
        set_det_coefficients(csfs, det_csf)

    # reading orbitals coefficients
    while 'ORBITALS IN PRIMITIVE BASIS FUNCTIONS' not in line:
        line = xmo_file.readline()
        if not line:
            sys.exit("Error: orbitals in primitive basis functions not found in .xmo file")
    for i in range(4 + 2 * xmvb3):
        line = xmo_file.readline()

    number_orbitals = len(orbitals)
    if wf_type == 'orb':
        for i in range(number_basisfunctions - number_orbitals):
            orbitals.append(Orbital())

    for i in range((len(orbitals) + 4) // 5):
        orbitals_to_read = min(len(orbitals) - i * 5, 5)
        for j in range(number_basisfunctions):
            words = line.split()
            for k in range(orbitals_to_read):
                orbitals[5 * i + k].coefficients.append(float(words[k - orbitals_to_read]))
            line = xmo_file.readline()
        for j in range(2 + 2 * xmvb3):
            line = xmo_file.readline()

    xmo_file.close()

    wf = WaveFunction(input_name, orbital_format, basis, charge, multiplicity, atoms, orbitals,
                      csfs, Jastrow(), total_symmetry, fragments)

    return wf


def set_det_coefficients(csfs, det_csf):
    csf_list = list(range(len(csfs)))
    while len(csf_list) != 0:
        unique_det_in_csf = False
        for csf_index in csf_list:
            this_csf = csfs[csf_index]

            # determining unique determinant
            this_det = Determinant()
            for det in this_csf.determinants:
                for k in csf_list:
                    if k != csf_index and det in csfs[k].determinants:
                        break
                else:  # no break
                    this_det = det
                    unique_det_in_csf = True
                    break

            if unique_det_in_csf:
                # initializing the det coefficient with the read det coefficient
                # index() works, since determinants are defined as __eq__ regardless of the coefficient
                det_coefficient = det_csf.determinants[
                    det_csf.determinants.index(this_det)].coefficient

                # substracting coefficient contributions from csfs that are already out of csf_list
                for other_csf_index, other_csf in enumerate(csfs):
                    if not (other_csf_index in csf_list or other_csf_index == csf_index):
                        if this_det in other_csf.determinants:
                            other_det_index = other_csf.determinants.index(this_det)
                            det_coefficient -= other_csf.determinants[other_det_index].coefficient \
                                               * other_csf.coefficient

                # setting all det coefficients of this_csf
                for det in this_csf.determinants:
                    det.coefficient /= abs(det.coefficient)
                    det.coefficient *= det_coefficient / csfs[csf_index].coefficient

                # removing csf_index from csf_list
                del csf_list[csf_list.index(csf_index)]

                # start screening csf_list from 0th entry again
                break

        else:  # no break
            sys.exit('Error: This should not happen, no csf with unique det found.')


def read_atoms(bohr, charge, info_file, line, number_atoms):
    atoms = []
    while '.0' not in line:
        line = info_file.readline()
        if not line:
            sys.exit('Error: no geometry found in INFO file.')
    for i in range(number_atoms):
        atom = Atom()
        words = line.split()
        atom.atomic_number = [item.upper() for item in periodic_table].index(words[0])
        for j in range(3):
            atom.position.append(float(words[j + 2].replace("D", "E")) * bohr)
        charge += int(float(words[1].replace("D", "E")))
        atoms.append(atom)
        line = info_file.readline()
    return atoms, charge, line
