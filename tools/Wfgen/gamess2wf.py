#!/usr/bin/env python

"""
 Copyright (C) 2022 Leonard Reuter

SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
if sys.version_info[0] < 3:
    sys.exit('This script requires Python 3')

from .Utils import *


def gamess_in(input_name,dat_name,basis,wf_type):
    out_filename = input_name
    dat_filename = dat_name
    orbital_format = 'gms'

    #opening files
    try:
        out_file = open(out_filename)
    except FileNotFoundError:
        not_found(out_filename)

    try:
        dat_file = open(dat_filename)
    except FileNotFoundError:
        not_found(dat_filename)

    line = out_file.readline()
    while "THE POINT GROUP OF THE MOLECULE IS" not in line:
        line = out_file.readline()
    total_symmetry = line.split()[-1]
    if 'N' in total_symmetry:
        while "THE ORDER OF THE PRINCIPAL AXIS IS" not in line:
            line = out_file.readline()
        total_symmetry = total_symmetry.replace('N', line.split()[-1])
    symmetry_list = symmetries[total_symmetry]

    atoms = []
    while "ATOM      ATOMIC" not in line:
        line = out_file.readline()
    for _ in range(2):
        line = out_file.readline()
    words = line.split()
    upper_pse = [item.upper() for item in periodic_table]
    while len(words) > 1:
        atom = Atom()
        atom.atomic_number = upper_pse.index(words[0].upper())
        for j in range(3):
            atom.position.append(float(words[j+2]))
        atoms.append(atom)
        line = out_file.readline()
        words = line.split()

    #for atom in atoms:
    #    atom.write(sys.stdout, False)

    while "NUMBER OF ELECTRONS" not in line:
        line = out_file.readline()
    number_electrons = int(line.split()[-1])

    while "CHARGE OF MOLECULE" not in line:
        line = out_file.readline()
    charge = int(line.split()[-1])

    while "SPIN MULTIPLICITY" not in line:
        line = out_file.readline()
    multiplicity = int(line.split()[-1])

    # construction for wf_type = 'sd'
    csfs = []
    if wf_type == 'csf':
        while "DETERMINANT CONTRIBUTION TO CSF'S" not in line:
            line = out_file.readline()
            if not line:
                sys.exit('Error: no CSFs found in .out file')
        for _ in range(2):
            line = out_file.readline()
        while "CASE VECTOR" in line:
            csf = Csf()
            for _ in range(4):
                line = out_file.readline()
            while 'C(' in line:
                determinant = Determinant()
                words = line.split("=")[-1].split()
                determinant.coefficient = float(words[0])
                for _ in range(2):
                    del(words[0])

                for i in range((number_electrons - len(words))//2):
                    # adding core electrons
                    determinant.orbital_list += [i+1, -(i+1)]
                determinant.orbital_list += [int(word) for word in words]
                determinant.sort()
                csf.determinants.append(determinant)
                line = out_file.readline()
            csfs.append(csf)
        while "CSF      COEF" not in line:
            line = out_file.readline()
        for _ in range(2):
            line = out_file.readline()
        for i in range(len(csfs)):
            words = line.split()
            if i+1 != int(words[0]):
                print(f'Warning: coeff of CSF {i+1} is below gamess printtol and set to 0.0')
                csfs[i].coefficient = 0.0
            else:
                csfs[i].coefficient = float(words[1])
                csfs[i].occupation = list(words[2])
                line = out_file.readline()

    elif wf_type == 'sd':
        determinant = Determinant()
        for i in range((number_electrons + multiplicity - 1) // 2):
            determinant.orbital_list.append(i + 1)
        for i in range((number_electrons - multiplicity + 1) // 2):
            determinant.orbital_list.append(-(i + 1))
        csf = Csf()
        csf.determinants.append(determinant)
        csf.number_determinants = 1
        csfs.append(csf)

    out_file.close()

    line = dat_file.readline()
    while '$VEC' not in line:
        line = dat_file.readline()
        if not line:
            sys.exit('Error: no $VEC section found in .dat file')
    print('Warning: this script always uses the first $VEC section in the .dat file')
    line = dat_file.readline()
    lines = 0
    number_orbitals = 0
    while '$END' not in line:
        number_orbitals = int(line.split()[0])
        line = dat_file.readline()
        lines += 1

    dat_file.seek(0)
    while '$VEC' not in line:
        line = dat_file.readline()
    line = dat_file.readline()
    mo_lines = lines//number_orbitals

    orbitals = []
    for i in range(number_orbitals):
        orbital = Orbital()
        for j in range(mo_lines):
            for k in range(len(line)//15):
                orbital_string = ''
                for l in range(15):
                    orbital_string += line[k * 15 + l + 5]
                orbital.coefficients.append(float(orbital_string.replace("D", "E")))
            line = dat_file.readline()
        orbitals.append(orbital)
    
    number_basisfunctions = len(orbitals[0].coefficients)

    #for i, orbital in enumerate(orbitals):
    #    orbital.write(sys.stdout, i, 'gms')

    dat_file.close()

    jastrow = Jastrow()
    jastrow.type = 'none'

    print('Warning: Reading of orbital symmetries is not yet implemented for gamess')
    wf = WaveFunction(input_name.split('.')[0], orbital_format, basis, charge, multiplicity, atoms, orbitals,
                      csfs, jastrow, total_symmetry, symmetry_list, bohr=True)

    return wf
