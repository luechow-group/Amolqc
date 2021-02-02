#
""""
 Copyright (C) 2011 Arne Luechow

 SPDX-License-Identifier: GPL-3.0-or-later
"""
# pse module
#


psel = ['X' ,'H' ,'He',
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

def elem(idx):
    assert idx>0 and idx<len(psel),"pse.getElement: internal pse too small"
    return psel[idx]

def Z(elem):
    assert elem in psel,"pse.getZ: element is not in internal pse"
    return psel.index(elem)


if __name__ == "__main__":
    print("testing the pse module")
    assert psel[8] == 'O', "Test 1 failed"
    print(Z('Cl'))

    
