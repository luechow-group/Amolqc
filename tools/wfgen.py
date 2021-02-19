#!/usr/bin/env python

""""
 Copyright (C) 2017-2019 Leonard Reuter

 SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
from Wfgen.xmvb2wf import xmvb_in
from Wfgen.mol2wf import molpro_in
from Wfgen.orca2wf import orca_in
from Wfgen.amolqc2wf import amolqc_in

if len(sys.argv) == 1:
    sys.exit('''
    wfgen.py reads the wave function from the output files of XMVB or Molpro
    or an amolqc .wf file. It can write .wf files, punch MOs in .vec files,
    sort determinants, combine csfs ...
    
    usage:   |             read section               ||   command section   |
    wfgen.py program file1 (file2) (wave function type) command (command2) ...
    
    The wftype is 'det','csf' or 'sd'.
    'sd' (single determinant) only for Molpro.
    
    Detailled read input schemes:
        Molpro: wfgen.py molpro file.out file.molden basis wftype
        orca:   wfgen.py orca file.mkl file.out basis wftype
                note: generate mkl file from gbw file with orca_2mkl
        XMVB:   wfgen.py xmvb   file.xmo INFO        basis wftype
        Amolqc: wfgen.py amolqc file.wf
    
    Path to the files can be given as well.
    
    Commands are: (several commands can be used at once)
        write              - writes a .wf file
        punch_orbs         - punches orbitals in a .vec file
        sort_dets          - sorts the orbital list of every determinant
        sort_ci            - sorts the CI vector by coefficient
        convert            - converts csf type to det type
        cut_ci [threshold] - cuts csf CI vector at given threshold
        add_symmetry       - adds symmetry (eg. delta for C2v)
        symm_combine       - combines csfs based on symmetry (molpro only)
        coeff_combine      - combines csfs with equal (diff < 1E-6) coeffs
        rename [name]      - .wf and .vec files written will have this name
        alt_sort_ci        - sorts the CI vector by occupation
        check_norm         - checks normalization
        normalize          - normalizes the determinant part of the wave function
        cut_orbs           - removes virtual orbitals (not in the active space)
        excitation         - calculates the sum of the squares of coefficients
                             of determinants/csfs, where given orbitals are not
                             doubly occupied
        punch_moopt        - punches the orbital optimization section for Amolqc input
                             and deletes all orbitals, that are not required
        max_dets           - prints the maximum number of determinants per csf
        share_dets         - prints indices of csfs that 'share' determinants
        replace_orbs       - replaces orbital indices in determinants
        reorder_atoms      - reorders atoms (and orbitals accordingly)
        separate_electrons - separates electrons (e.g. for pi-only calculations)
        
    ''')

program = sys.argv[1]
arguments = 0

if program == 'xmvb':
    arguments = 6
    input_name = sys.argv[2].split('.xmo')[0]
    info_name = sys.argv[3]
    basis = sys.argv[4]
    wf_type = sys.argv[5]
    if not (wf_type in ['csf', 'det', 'orb']):
        sys.exit("Error: wave function type '"+wf_type+"' is not supported!")
    
    wf = xmvb_in(input_name, info_name, basis, wf_type)
    if wf.get_csfs_sharing_determinants(verbose=False):
        print('ATTENTION: csfs sharing determinants detected, wave function may be wrong')

elif program == 'molpro':
    arguments = 6
    input_name = sys.argv[2].split('.out')[0]
    molden_name = sys.argv[3].split('.molden')[0]
    basis = sys.argv[4]
    wf_type = sys.argv[5]
    if not (wf_type == 'csf' or wf_type == 'det' or wf_type == 'sd'):
        sys.exit("Error: wave function type '" + wf_type + "' is unknown!")
    wf = molpro_in(input_name, molden_name, basis, wf_type)
    
elif program == 'orca':
    arguments = 6
    mkl_name = sys.argv[2]
    orcaout_name = sys.argv[3]
    basis = sys.argv[4]
    wf_type = sys.argv[5]
    if wf_type not in ['sd', 'csf']:
        sys.exit("Error: wave function type '" + wf_type + "' is unknown!")
    wf = orca_in(mkl_name, orcaout_name, basis, wf_type)
    
elif program == 'amolqc':
    arguments = 3
    input_name = sys.argv[2].split('.wf')[0]
    wf = amolqc_in(input_name)

else:
    sys.exit("Error: program '"+program+"' unknown.")

i = 0

written = False
while i < len(sys.argv) - arguments:
    command = sys.argv[i + arguments]
    if command == 'write':
        if program == 'amolqc' and 'rename' not in sys.argv:
            wf.title += '_new'
            print('_new has been appended to the .wf filename in order to avoid overwriting.')
        if wf.charge != 0 and not wf.atomic_charges():
            wf.set_charges()
        wf.write()
        written = True
    elif command == 'punch_orbs':
        wf.punch_vec()
    elif command == 'sort_dets':
        wf.sort_dets()
    elif command == 'convert':
        wf.convert_to_det()
    elif command == 'sort_ci':
        wf.sort_ci()
    elif command == 'cut_ci':
        i += 1
        threshold = float(sys.argv[i + arguments])
        wf.cut_ci(threshold)
    elif command == 'rename':
        i += 1
        name = sys.argv[i + arguments]
        wf.rename(name)
    elif command == 'alt_sort_ci':
        wf.alt_sort_ci()
    elif command == 'check_norm':
        wf.check_normalization()
    elif command == 'excitation':
        wf.check_excitation()
    elif command == 'normalize':
        wf.normalize()
    elif command == 'coeff_combine':
        wf.coeff_combine()
    elif command == 'symm_combine':
        if not wf.symmetry:
            sys.exit('Error: symm_combine called, but no symmetry found. Has to be used directly when reading in molpro wf.')
        if wf.symmetry == 'VB':
            sys.exit('Error: symm_combine is not yet implemented for VB wave functions.')
        wf.symm_combine()
    elif command == 'cut_orbs':
        wf.cut_orbs()
    elif command == 'max_dets':
        wf.max_dets()
    #elif command == 'force_combine':
    #    wf.force_combine()
    elif command == 'add_symmetry':
        wf.add_symmetry()
    elif command == 'punch_moopt':
        wf.punch_moopt()
    elif command == 'share_dets':
        wf.get_csfs_sharing_determinants()
    elif command == 'replace_orbs':
        wf.replace_orbital_indices()
    elif command == 'reorder_atoms':
        wf.reorder_atoms()
    elif command == 'separate_electrons':
        wf.separate_electrons()
    elif command == 'print_csf':
        i += 1
        index = int(sys.argv[i + arguments])
        wf.print_csf(index)
    else:
        sys.exit('Error: command '+command+' unknown!')
    i += 1

if not written:
    print("Warning: no wave function file has been written, since the command 'write' was not used")
