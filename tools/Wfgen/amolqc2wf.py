#!/usr/bin/env python

"""
 Copyright (C) 2017-2018 Leonard Reuter

SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
from .Utils import *


def read_dets(number_determinants,number_electrons,number_alpha_electrons,wf_file):
    determinants = []
    for i in range(number_determinants):
        line = wf_file.readline()
        words = line.split()
        determinant = Determinant()
        determinant.coefficient = float(words[0])
        for j in range(number_alpha_electrons):
            determinant.orbital_list.append(int(words[j+1]))
        for j in range(number_electrons - number_alpha_electrons):
            determinant.orbital_list.append(-int(words[j+1+number_alpha_electrons]))
        determinants.append(determinant)
    return determinants


def amolqc_in(input_name):
    wf_filename = input_name + '.wf'
    #opening file
    try:
        wf_file = open(wf_filename)
    except FileNotFoundError:
        not_found(wf_filename)

    #reading general section
    for i in range(2):
        line = wf_file.readline()
    general = {}
    general['atomic_charges'] = False
    general['evfmt'] = 'gau'
    general['geom'] = 'ang'
    while "$end" not in line:
        while line[0] == ' ':
            line = line[1:]
        words = line[:len(line)-1].split(', ')
        for i in range(len(words)):
            if '=' in words[i]:
                general[words[i].split('=')[0]] = words[i].split('=')[1]
            else:
                general[words[i]] = True
        line = wf_file.readline()

    charge = int(general['charge'])
    multiplicity = int(general['spin'])
    atomic_charges = general['atomic_charges']
    basis = general['basis']
    orbital_format = general['evfmt']
    title = general['title']
    jastrow_type = general['jastrow']
    bohr = False
    if general['geom'] == 'bohr':
        bohr = True
    
    #reading geometry
    while '$geom' not in line:
        line = wf_file.readline()
        if not line:
            sys.exit('Error: geom section not found')
    line = wf_file.readline()
    number_atoms = int(line)
    atoms = []
    for i in range(number_atoms):
        line = wf_file.readline()
        words = line.split()
        atom = Atom()
        atom.atomic_number = periodic_table.index(words[0])
        for j in range(3):
            atom.position.append(float(words[j+1]))
        if atomic_charges:
            atom.charge = int(words[4])
        atoms.append(atom)

    #reading jastrow
    jastrow = Jastrow()
    jastrow.type = jastrow_type
    if jastrow.type != 'none':
        while '$jastrow' not in line:
            line = wf_file.readline()
            if not line:
                sys.exit('Error: jastrow section not found')
        line = wf_file.readline()
        while '$end' not in line:
            jastrow.lines.append(line)
            line = wf_file.readline()

    #reading orbitals
    while '$mos' not in line:
        line = wf_file.readline()
        if not line:
            sys.exit('Error: mo section not found')

    line = wf_file.readline()
    words = line.split()
    number_orbitals = int(words[0])
    for i in range(2):
        line = wf_file.readline()
    lines = 0
    while '$end' not in line:
        line = wf_file.readline()
        lines += 1
    wf_file.seek(0)
    mo_lines = lines//number_orbitals
    while '$mos' not in line:
        line = wf_file.readline()
    for i in range(3):
        line = wf_file.readline()
    orbitals = []
    if orbital_format in ['gau', 'fre']:
        mo_lines -= 1

    for i in range(number_orbitals):
        orbital = Orbital()
        if orbital_format in ['gau', 'fre']:
            line = wf_file.readline()
        if orbital_format == 'fre':
            for j in range(mo_lines):
                for k in range(len(line.split())):
                    orbital_string = line.split()[k]
                    orbital.coefficients.append(float(orbital_string.replace("e", "E")))
                line = wf_file.readline()
        elif orbital_format in ['gau', 'gms']:
            for j in range(mo_lines):
                for k in range(len(line)//15):
                    orbital_string = ''
                    for l in range(15):
                        if orbital_format == 'gau':
                            orbital_string += line[k * 15 + l]
                        elif orbital_format == 'gms':
                            orbital_string += line[k * 15 + l + 5]
                    orbital.coefficients.append(float(orbital_string.replace("D", "E")))
                line = wf_file.readline()
        else:
            sys.exit('Error: orbital format mol is not implemented for mo read')
        orbitals.append(orbital)

    number_basisfunctions = len(orbitals[0].coefficients)

    csfs = []
    
    #reading determinants
    dets_section = True
    while '$dets' not in line:
        line = wf_file.readline()
        if not line:
            dets_section = False
            wf_file.seek(0)
            break
    if dets_section:
        csf = Csf()
        line = wf_file.readline()
        if 'single' not in line:
            number_determinants = int(line)
            pos = wf_file.tell()
            line = wf_file.readline()
            words = line.split()
            number_electrons = len(words) - 1
            number_alpha_electrons = (multiplicity + number_electrons - 1)//2
            wf_file.seek(pos)
            csf.determinants = read_dets(number_determinants,number_electrons,number_alpha_electrons,wf_file)
            csfs.append(csf)
        else:
            number_electrons = int(input('What is the total number of electrons? '))
            determinant = Determinant()
            for i in range((number_electrons + multiplicity - 1) // 2):
                determinant.orbital_list.append(i + 1)
            for i in range((number_electrons - multiplicity + 1) // 2):
                determinant.orbital_list.append(-(i + 1))
            csf = Csf()
            csf.determinants.append(determinant)
            csf.number_determinants = 1
            csfs.append(csf)

    #reading csfs
    else:
        while '$csfs' not in line:
            line = wf_file.readline()
            if not line:
                sys.exit('Error: neither csfs nor dets found')
        line = wf_file.readline()
        pos = wf_file.tell()
        if 'single' not in line:
            number_csfs = int(line)
            for i in range(2):
                line = wf_file.readline()
            words = line.split()
            number_electrons = len(words) - 1
            number_alpha_electrons = (multiplicity + number_electrons - 1)//2 
            wf_file.seek(pos)
            csfs = []
            for i in range(number_csfs):
                line = wf_file.readline()
                words = line.split()
                csf = Csf()
                csf.coefficient = float(words[0])
                number_determinants = int(words[1])
                csf.determinants = read_dets(number_determinants,number_electrons,number_alpha_electrons,wf_file)
                csfs.append(csf)
        else:
            number_electrons = int(input('What is the total number of electrons? '))
            determinant = Determinant()
            for i in range((number_electrons + multiplicity - 1) // 2):
                determinant.orbital_list.append(i + 1)
            for i in range((number_electrons - multiplicity + 1) // 2):
                determinant.orbital_list.append(-(i + 1))
            csf = Csf()
            csf.determinants.append(determinant)
            csf.number_determinants = 1
            csfs.append(csf)

    wf = WaveFunction(input_name, orbital_format, basis, charge, multiplicity, atoms,
                      orbitals, csfs, jastrow, bohr=bohr)
    
    return wf
