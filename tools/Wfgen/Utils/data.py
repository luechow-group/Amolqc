""""
 Copyright (C) 2017-2019 Leonard Reuter

SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
if sys.version_info[0] < 3:
    sys.exit('This script requires Python 3')

periodic_table = ['X' ,'H' ,'He',
                  'Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne',
                  'Na','Mg','Al','Si','P' ,'S' ,'Cl','Ar',
                  'K' ,'Ca','Sc','Ti','V' ,'Cr','Mn','Fe',
                  'Co','Ni','Cu','Zn','Ga','Ge','As','Se',
                  'Br','Kr','Rb','Sr','Y' ,'Zr','Nb','Mo',
                  'Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
                  'Sb','Te','I' ,'Xe','Cs','Ba','La','Ce',
                  'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',
                  'Ho','Er','Tm','Yb','Lu','Hf','Ta','W',
                  'Re','Os','Ir','Pt','Au','Hg','Tl','Pb',
                  'Bi','Po','At','Rn','Fr','Ra','Ac','Th',
                  'Pa','U' ,'Np','Pu','Am','Cm','Bk','Cf',
                  'Es','Fm','Md','No','Lr','Rf','Db','Sg',
                  'Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl',
                  'Mc','Lv','Ts','Og']

# multiplication table for ANY! molpro symmetry
multiplication_table = [
    [1,2,3,4,5,6,7,8],
    [2,1,4,3,6,5,8,7],
    [3,4,1,2,7,8,5,6],
    [4,3,2,1,8,7,6,5],
    [5,6,7,8,1,2,3,4],
    [6,5,8,7,2,1,4,3],
    [7,8,5,6,3,4,1,2],
    [8,7,6,5,4,3,2,1]]

# symmetries
symmetries = {
    'D2h':['Ag','B3u','B2u','B1g','B1u','B2g','B3g','Au'],
    'C2v':['A1','B1','B2','A2'],
    'C2h':['Ag','Au','Bu','Bg'],
    'D2':['A','B3','B2','B1'],
    'Cs':["A'","A''"],
    'C2':['A','B'],
    'Ci':['Ag','Au'],
    'C1':['A']}
