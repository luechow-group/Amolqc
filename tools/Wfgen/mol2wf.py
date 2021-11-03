#!/usr/bin/env python

"""
 Copyright (C) 2017-2019 Leonard Reuter

SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
from fractions import Fraction
from .Utils import *


def molpro_in(input_name,molden_name,basis,wf_type):
    out_filename = input_name + '.out'
    molden_filename = molden_name+'.molden'
    orbital_format = 'gau'

    #opening files
    try:
        out_file = open(out_filename)
    except FileNotFoundError:
        not_found(out_filename)

    try:
        out_file.readlines()
    except UnicodeDecodeError:
        print('Warning: There is something weird with the molpro output:')
        out_file.close()
        out_file = open(out_filename, 'rb')
        lines = out_file.readlines()
        out_file.close()

        out_filename = out_filename.split('.out')[0]+'_new.out'
        out_file = open(out_filename, 'wb')
        for i, byte_line in enumerate(lines):
            line_is_fine = True
            try:
                byte_line.decode('ASCII')
            except UnicodeDecodeError:
                line_is_fine = False
                print(f'line {i} is broken')
            if line_is_fine:
                out_file.write(byte_line)

        print(f'Wrote a fixed .out file to {out_filename}')
        del lines

    out_file.close()
    out_file = open(out_filename, 'r')

    try:
        molden_file = open(molden_filename)
    except FileNotFoundError:
        not_found(molden_filename)

    #checking if basf order in molden is according to atom indices
    line = molden_file.readline()
    while '[GTO]' not in line:
        line = molden_file.readline()
        if not line:
            sys.exit('Error: gtos not found in moldenfile')
    gto_order = []
    line = molden_file.readline()
    nuclear_charges = []
    i = 0
    gtos = []
    while True:
        gto_order.append(int(line.split()[0]))
        nuclear_charges.append(line.split()[1])
        gtos.append('')
        while True:
            line = molden_file.readline()
            gtos[i] += line
            if len(line.split()) == 0:
                line = molden_file.readline()
                i += 1
                break
        if len(line.split()) == 0:
            break
    correct_order = True
    molden_file.seek(0)
    for i in range(len(gto_order)):
        if gto_order[i] != i + 1:
            correct_order = False
            break
    if not correct_order:
        molden_file.close()
        molden_file = open(molden_name+'_new.molden', 'w')
        old_molden_file = open(molden_filename)
        line = ''
        while '[Atoms]' not in line:
            line = old_molden_file.readline()
            molden_file.write(line)
        line = old_molden_file.readline()
        words = line.split()
        unsorted_atoms = []
        while len(words) == 6:
            atom = Atom()
            atom.atomic_number = [item.upper() for item in periodic_table].index(words[0].upper())
            atom.core_charge = int(words[2])
            for i in range(3):
                atom.position.append(float(words[3 + i]))
            unsorted_atoms.append(atom)
            line = old_molden_file.readline()
            words = line.split()
        sorted_atoms = []
        for i in range(len(gto_order)):
            sorted_atoms.append(unsorted_atoms[gto_order[i] - 1])
        for i in range(len(sorted_atoms)):
            sorted_atoms[i].write_molpro(molden_file,i)
        molden_file.write('[GTO]'+'\n')
        for i in range(len(gtos)):
            molden_file.write(str(i + 1).rjust(4)
                              + ' '
                              + nuclear_charges[i] + '\n'
                              + gtos[i])
        molden_file.write(' ' + '\n')

        while '[MO]' not in line:
            line = old_molden_file.readline()
        while line:
            molden_file.write(line)
            line = old_molden_file.readline()

        old_molden_file.close()
        molden_file.close()
        molden_file = open(molden_name+'_new.molden')
        print(f"Order of atoms in molden file edited and/or dummy atoms removed: {molden_name}_new.molden. ")

    #reading data from molden file
    line = molden_file.readline()

    while 'NELEC=' not in line:
        line = molden_file.readline()
        if not line:
            sys.exit('Error: number of electrons not found in moldenfile')
    words = line.split()
    number_electrons = int(float(words[1]))

    while '!SPIN=' not in line:
        line = molden_file.readline()
        if not line:
            sys.exit('Error: multiplicity not found in moldenfile')
    words = line.split()
    multiplicity = int(float(words[1])) + 1

    charge = -1 * number_electrons
    number_atoms = 0
    atoms = []

    while '[Atoms]' not in line:
        line = molden_file.readline()
        if not line:
            sys.exit('Error: geometry not found in moldenfile')
    line = molden_file.readline()
    words = line.split()
    while len(words) == 6:
        charge += int(words[2])
        number_atoms += 1
        atom = Atom()
        atom.atomic_number = [item.upper() for item in periodic_table].index(words[0].upper())
        for i in range(3):
            atom.position.append(float(words[3 + i]))
        atoms.append(atom)
        line = molden_file.readline()
        words = line.split()

    #determining number of basis functions

    while '[MO]' not in line:
        line = molden_file.readline()
        if not line:
            sys.exit('Error: molecular orbitals not found in moldenfile')
    for i in range(5):
        line = molden_file.readline()
    words = line.split()
    number_basisfunctions = 0
    while words[0] != 'Sym=':
        number_basisfunctions += 1
        line = molden_file.readline()
        words = line.split()

    #reading MOs
    number_orbitals = 0
    orbitals = []
    orbital_symmetries = []

    molden_file.seek(0)
    while '[MO]' not in line:
        line = molden_file.readline()
    line = molden_file.readline()
    while line:
        orbital = Orbital()
        number_orbitals += 1
        words = line.split()
        orbital_symmetries.append(words[1])
        orbital.symmetry = words[1]
        for i in range(4):
            line = molden_file.readline()
        for i in range(number_basisfunctions):
            words = line.split()
            orbital.coefficients.append(float(words[1]))
            line = molden_file.readline()
        orbitals.append(orbital)
    molden_file.close()
    
    csfs = []
    total_symmetry = None
    symmetry_list = []

    line = out_file.readline()
    while 'Point group' not in line:
        line = out_file.readline()
        if not line:
            sys.exit('Error: point group could not be read!')
    total_symmetry = line.split()[2]
    symmetry_list = symmetries[total_symmetry]

    if wf_type != 'sd':
        #reading the symmetry of the occupied orbitals
        number_irreps = len(symmetry_list)

        #reading determinants
        if wf_type == 'det':
            sys.exit("wf_type 'det' is currently not supported")

        #reading csfs
        if wf_type == 'csf':
            while 'Reference space' not in line:
                line = out_file.readline()
                if not line:
                    sys.exit('Error: csfs of mrci not found!')

            while 'orbitals:' not in line:
                line = out_file.readline()

            # read number and symmetry of core orbitals
            number_core_orbitals = 0;
            core_orbital_symmetries = [0]*number_irreps
            if 'core' in line:
                words = line.split()
                number_core_orbitals = int(words[4])
                for i in range(6):
                    del words[0]

                del words[-1]
                if len(words) != number_irreps:
                    while len(words) < number_irreps:
                        words.append('0')

                for i in range(len(words)):
                    core_orbital_symmetries[i] = int(words[i])
                line = out_file.readline()

            # read number and symmetry of closed orbitals
            number_closed_orbitals = 0;
            closed_orbital_symmetries = [0]*number_irreps
            if 'closed' in line:
                words = line.split()
                number_closed_orbitals = int(words[4])
                for i in range(6):
                    del words[0]

                del words[-1]
                if len(words) != number_irreps:
                    while len(words) < number_irreps:
                        words.append('0')

                for i in range(len(words)):
                    closed_orbital_symmetries[i] = int(words[i])
                line = out_file.readline()

            # read number and symmetry of active orbitals
            number_active_orbitals = 0
            active_orbital_symmetries = [0]*number_irreps
            if 'active' in line:
                words = line.split()
                number_active_orbitals = int(words[4])
                for i in range(6):
                    del words[0]

                del words[-1]
                if len(words) != number_irreps:
                    while len(words) < number_irreps:
                        words.append('0')

                for i in range(len(words)):
                    active_orbital_symmetries[i] = int(words[i])
                line = out_file.readline()

            # calculate number of orbitals and electrons and symmetries of occupation vector string
            number_vector_orbitals = number_active_orbitals + number_closed_orbitals
            number_vector_electrons = number_electrons - 2 * number_core_orbitals
            vector_orbital_symmetries = []
            for i in range(number_irreps):
                vector_orbital_symmetries.append(closed_orbital_symmetries[i]
                                                 + active_orbital_symmetries[i])

            # build orbital map with 'orbital_symmetries' from .molden file
            orbital_map = []
            orbital_symmetry_counter = core_orbital_symmetries
            for i in range(number_irreps):
                for j in range(vector_orbital_symmetries[i]):
                    orbital_symmetry_counter[i] += 1
                    symmetry_string = str(orbital_symmetry_counter[i]) + '.' + str(i + 1)
                    index = orbital_symmetries.index(symmetry_string)
                    orbital_map.append(index + 1)

            # build csfs
            csf_coefficients = []
            csf_occupations = []
            occupation_scheme = ''
            spin_function = []
            spin_coefficient = ''

            # for speedup: determine number of maximum singly occupied orbitals
            max_singly_occ = min(2 * number_vector_orbitals - number_vector_electrons, number_vector_electrons)

            # build genealogical_spin_functions (call genealogical.py)
            spin = Fraction(multiplicity - 1, 2)
            if spin == 0:
                X = create_genealogical_spin_functions(max_singly_occ, singlet=True)
            else:
                X = create_genealogical_spin_functions(max_singly_occ)

            # read csf coefficients and occupation vector strings
            while 'Reference coefficients' not in line:
                line = out_file.readline()
                if not line:
                    sys.exit('Error: csfs not found!')
            for i in range(2):
                line = out_file.readline()
            words = line.split()

            number_of_states = len(words) - 1
            if number_of_states > 1:
                print('state-averaged wave function detected')
                state = int(input('give the number of the state you want to read: '))
            else:
                state = 1

            while words:
                csf_coefficients.append(float(words[state]))
                csf_occupations.append(list(words[0]))
                line = out_file.readline()
                words = line.split()

            # build csfs with spin functions
            for i in range(len(csf_coefficients)):
                csf = Csf()
                csf.occupation = ['0'] * number_orbitals

                # put 2 electrons in all core orbitals
                for j in range(number_core_orbitals):
                    csf.occupation[j] = '2'
                # add strings from the csf occupation vector string
                for j in range((len(orbital_map))):
                    csf.occupation[orbital_map[j] - 1] = csf_occupations[i][j]
                # add coefficient
                csf.coefficient = csf_coefficients[i]

                # find occupation scheme (only singly occupied orbitals)
                for j in range(number_vector_orbitals):
                    if csf_occupations[i][j] == '/':
                        occupation_scheme += '/'
                    elif csf_occupations[i][j] == '\\':
                        occupation_scheme += '\\'

                if occupation_scheme != '':  # if singly occupied orbitals in csf
                    Xlist = X[len(occupation_scheme)][spin//1][int(2*spin)]
                    for sf in Xlist:
                        words = str(sf).split()
                        if words[0] == occupation_scheme:
                            spin_function = list(str(sf).split())
                    del spin_function[0]

                    for j in range(len(spin_function)):
                        occupation = list(csf_occupations[i])
                        single_occupation = []
                        letters = list(spin_function[j])
                        for k in range(len(spin_function[j])):
                            if letters[k] in ['a', 'b']:
                                single_occupation.append(letters[k])
                            elif letters[k] != '*':
                                spin_coefficient += letters[k]
                        coefficient = float(abs(Fraction(spin_coefficient)))**0.5 \
                                      * int(Fraction(spin_coefficient)/abs(Fraction(spin_coefficient)))
                        for k in range(len(occupation)):
                            if occupation[k] in ['/', '\\']:
                                occupation[k] = single_occupation[0]
                                del single_occupation[0]
                        determinant = build_det(number_core_orbitals, occupation, orbital_map, coefficient)
                        csf.determinants.append(determinant)
                        spin_coefficient = ''


                else:  # if no singly occupied orbitals in csf
                    occupation = list(csf_occupations[i])
                    coefficient = 1.0
                    determinant = build_det(number_core_orbitals,occupation,orbital_map,coefficient)
                    csf.determinants.append(determinant)

                csfs.append(csf)
                occupation_scheme = ''

    # construction for wf_type = 'sd'
    else:
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
    jastrow = Jastrow()
    jastrow.type = 'none'

    wf = WaveFunction(input_name, orbital_format, basis, charge, multiplicity, atoms, orbitals,
                      csfs, jastrow, total_symmetry, symmetry_list)
    return wf
