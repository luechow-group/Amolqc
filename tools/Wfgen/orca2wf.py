""""
 Copyright (C) 2020 Arne Luechow

SPDX-License-Identifier: GPL-3.0-or-later
"""

import numpy as np
import sys

from .Utils import *
from fractions import Fraction
#
# transformation of orca mkl file to amolqc wf.
# use transformation dictionary as constructed from 'pure_and_cartesian_gaussians.py'
# using definition of real solid harmonics from Helgaker, or recursion formula
# agreement with GAMESS (ISPHER=1) only assuming that orca uses _opposite sign_ for
# f3c and f3s, g3c, g3s, g4c, g4s
# Arne Luechow, 2020

orca_gms_transformation_dict = {
    'd0c' : [-1/2, -1/2, 1, 0, 0, 0],
    'd1c' : [0, 0, 0, 0, 1, 0],
    'd1s' : [0, 0, 0, 0, 0, 1],
    'd2c' : [np.sqrt(3)/2, -np.sqrt(3)/2, 0, 0, 0, 0],
    'd2s' : [0, 0, 0, 1, 0, 0],
    'f0c' : [0, 0, 1, 0, -3*np.sqrt(5)/10, 0, -3*np.sqrt(5)/10, 0, 0, 0],
    'f1c' : [-np.sqrt(6)/4, 0, 0, 0, 0, -np.sqrt(30)/20, 0, np.sqrt(30)/5, 0, 0],
    'f1s' : [0, -np.sqrt(6)/4, 0, -np.sqrt(30)/20, 0, 0, 0, 0, np.sqrt(30)/5, 0],
    'f2c' : [0, 0, 0, 0, np.sqrt(3)/2, 0, -np.sqrt(3)/2, 0, 0, 0],
    'f2s' : [0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    'f3c' : [-np.sqrt(10)/4, 0, 0, 0, 0, 3*np.sqrt(2)/4, 0, 0, 0, 0],
    'f3s' : [0, np.sqrt(10)/4, 0, -3*np.sqrt(2)/4, 0, 0, 0, 0, 0, 0],
    'g0c' : [3/8, 3/8, 1, 0, 0, 0, 0, 0, 0, 3*np.sqrt(105)/140, -3*np.sqrt(105)/35, -3*np.sqrt(105)/35, 0, 0, 0],
    'g1c' : [0, 0, 0, 0, -3*np.sqrt(70)/28, 0, 0, np.sqrt(70)/7, 0, 0, 0, 0, 0, -3*np.sqrt(14)/28, 0],
    'g1s' : [0, 0, 0, 0, 0, 0, -3*np.sqrt(70)/28, 0, np.sqrt(70)/7, 0, 0, 0, -3*np.sqrt(14)/28, 0, 0],
    'g2c' : [-np.sqrt(5)/4, np.sqrt(5)/4, 0, 0, 0, 0, 0, 0, 0, 0, 3*np.sqrt(21)/14, -3*np.sqrt(21)/14, 0, 0, 0],
    'g2s' : [0, 0, 0, -np.sqrt(35)/14, 0, -np.sqrt(35)/14, 0, 0, 0, 0, 0, 0, 0, 0, 3*np.sqrt(7)/7],
    'g3c' : [0, 0, 0, 0, -np.sqrt(10)/4, 0, 0, 0, 0, 0, 0, 0, 0, 3*np.sqrt(2)/4, 0],
    'g3s' : [0, 0, 0, 0, 0, 0, np.sqrt(10)/4, 0, 0, 0, 0, 0, -3*np.sqrt(2)/4, 0, 0],
    'g4c' : [-np.sqrt(35)/8, -np.sqrt(35)/8, 0, 0, 0, 0, 0, 0, 0, 3*np.sqrt(3)/4, 0, 0, 0, 0, 0],
    'g4s' : [0, 0, 0, -np.sqrt(5)/2, 0, np.sqrt(5)/2, 0, 0, 0, 0, 0, 0, 0, 0, 0]    
}

def orca_in(mkl_name, out_name, basis, wf_type):
    mol = RichOrbitals()
    mol.read_orca_MKL_file(mkl_name)
    mol.add_cartesian_coefs(orca_gms_transformation_dict)
    input_name = mkl_name.split(".mkl")[0]
    orbital_format = 'gms'
    charge, mult = mol.charge_and_mult()
    number_electrons = mol.number_of_electrons()

    number_orbitals = mol.number_of_mos()
    orbitals = []
    orbital_symmetries = []
    for i in range(number_orbitals):
        orbital = Orbital()
        for a in mol.atoms:
            for bf in a.basis_functions:
                for ss in mol.gms_cartorder[bf.typ]:
                    orbital.coefficients.append(bf.mocoef[ss][i])
        orbitals.append(orbital)

    csfs = []
    total_symmetry = None
    symmetry_list = []

    with open(out_name, 'r') as outfile:
        line = outfile.readline()

        if wf_type != 'sd':
            #reading csfs
            if wf_type == 'csf':
                while 'Extended CI Printing' not in line:
                    line = outfile.readline()
                    if not line:
                        sys.exit('Error: Configurations not found')
                for _ in range(5):
                    line = outfile.readline()

                words = line.split()
                occupation = list(words[1])
                number_active_orbitals = len(occupation)
                number_active_electrons = sum([int(element) for element in occupation])
                number_core_electrons = number_electrons - number_active_electrons
                assert(number_core_electrons % 2 == 0), 'Number of core electrons has to be even'
                number_core_orbitals = number_core_electrons // 2
                orbital_map = list(range(number_core_orbitals+1,number_active_orbitals+number_core_orbitals+1))

                # for speedup: determine number of maximum singly occupied orbitals
                max_singly_occ = min(2 * number_active_orbitals - number_active_electrons, number_active_electrons)

                # build genealogical_spin_functions (call genealogical.py)
                spin = Fraction(mult - 1, 2)
                if spin == 0:
                    X = create_genealogical_spin_functions(max_singly_occ, singlet=True)
                else:
                    X = create_genealogical_spin_functions(max_singly_occ)

                while 'CFG' in line:
                    words = line.split()
                    cfg_occupation = list(words[1])
                    number_singly_occupied = cfg_occupation.count('1')
                    i = 1
                    for _ in range(3):
                        line = outfile.readline()
                    while 'CSF' in line:
                        csf = Csf()
                        words = line.split()
                        csf.coefficient = float(words[-1])

                        if number_singly_occupied != 0:
                            spin_function = list(str(X[number_singly_occupied][spin//1][int(2*spin)][-i]).split())
                            del spin_function[0]

                            for j in range(len(spin_function)):
                                occupation = list(cfg_occupation)
                                single_occupation = []
                                spin_coefficient = ''
                                letters = list(spin_function[j])
                                for k in range(len(spin_function[j])):
                                    if letters[k] in ['a', 'b']:
                                        single_occupation.append(letters[k])
                                    elif letters[k] != '*':
                                        spin_coefficient += letters[k]
                                coefficient = float(abs(Fraction(spin_coefficient)))**0.5 \
                                              * int(Fraction(spin_coefficient)/abs(Fraction(spin_coefficient)))
                                for k in range(len(occupation)):
                                    if occupation[k] == '1':
                                        occupation[k] = single_occupation[0]
                                        del single_occupation[0]
                                determinant = build_det(number_core_orbitals, occupation, orbital_map, coefficient)
                                csf.determinants.append(determinant)
                        else:
                            determinant = build_det(number_core_orbitals, cfg_occupation, orbital_map, 1.0)
                            csf.determinants.append(determinant)

                        i += 1
                        csfs.append(csf)
                        for _ in range(2):
                            line = outfile.readline()
        else:
            determinant = Determinant()
            for i in range((number_electrons + mult - 1) // 2):
                determinant.orbital_list.append(i + 1)
            for i in range((number_electrons - mult + 1) // 2):
                determinant.orbital_list.append(-(i + 1))
            csf = Csf()
            csf.determinants.append(determinant)
            csf.number_determinants = 1
            csfs.append(csf)

    outfile.close()
    jastrow = Jastrow()
    jastrow.type = 'none'

    wf = WaveFunction(input_name, orbital_format, basis, charge, mult, mol.atoms, orbitals,
                      csfs, jastrow, total_symmetry, symmetry_list)

    return wf










