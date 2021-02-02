#!/usr/bin/env python
""""
 Copyright (C) 2011 Arne Luechow

 SPDX-License-Identifier: GPL-3.0-or-later
"""
import sys
from Basis import basisset

if len(sys.argv)<5:
    print("""usage: runConvertBasisFile.py intype infile outtype outfile [cvt [key]]
    where intype:'gbs' or 'emsl' or 'abs' and outtype:'gbs', 'gms' or 'abs'
    gbs is the gaussian basis set format, gms the gamess format, abs the amolqc format
    infile and outfile refer to the filenames
    cvt is a sto2gto conversion file required for STO -> GTO conversion
    key is an optional key for selection of the GTO expansion""")
else:
    if len(sys.argv)>6:
        cvt = sys.argv[5]
        key = sys.argv[6]
    elif len(sys.argv)>5:
        cvt = sys.argv[5]
        key = None
    else:
        cvt = None
        key = None
    try:
        basisset.convertBasisFile(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],cvt,key)
    except ValueError as e:
        print("An error occurred:",e)
        
 
    
    

