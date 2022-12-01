#!/usr/bin/env python

""""
 Copyright (C) 2020 Arne Luechow

SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
if sys.version_info[0] < 3:
    sys.exit('This script requires Python 3')

#
import getopt
from .data import periodic_table
from .atom import Atom
from .basisFunction import BasisFunction


class RichOrbitals:
    """class containing all mo coeffs with atom and basis information"""

    def __init__(self):
        self.orca_basisorder = {
            's' : ['s'],
            'p' : ['x', 'y', 'z'],
            'd' : ['d0c', 'd1c', 'd1s', 'd2c', 'd2s'],
            'f' : ['f0c', 'f1c', 'f1s', 'f2c', 'f2s', 'f3c', 'f3s'],
            'g' : ['g0c', 'g1c', 'g1s', 'g2c', 'g2s', 'g3c', 'g3s', 'g4c', 'g4s']
        }
        self.gms_cartorder = {
            's' : ['s'],
            'p' : ['x', 'y', 'z'],
            'd' : ['xx', 'yy', 'zz', 'xy', 'xz', 'yz'],
            'f' : ['xxx', 'yyy', 'zzz', 'xxy', 'xxz', 'yyx', 'yyz', 'zzx', 'zzy', 'xyz'],
            'g' : ['xxxx', 'yyyy', 'zzzz', 'xxxy', 'xxxz', 'yyyx', 'yyyz', 'zzzx', 'zzzy', 'xxyy', 'xxzz', 'yyzz', 'xxyz', 'yyxz', 'zzxy']
        }
        # molden order TODO: check!
        self.molden_cartorder = {
            's' : ['s'],
            'p' : ['x', 'y', 'z'],
            'd' : ['xx', 'yy', 'zz', 'xy', 'xz', 'yz'],
            'f' : ['xxx', 'yyy', 'zzz', 'xyy', 'xxy', 'xxz', 'xzz', 'yzz', 'yyz', 'xyz'],
            'g' : ['xxxx', 'yyyy', 'zzzz', 'xxxy', 'xxxz', 'yyyx', 'yyyz', 'zzzx', 'zzzy', 'xxyy', 'xxzz', 'yyzz', 'xxyz', 'yyxz', 'zzxy']
        }
        self.charge = 0
        self.mult = 1
        self.atoms = []

    def read_orca_MKL_file(self, filename):

        # read CHAR_MULT block
        try:
            file = open(filename, 'r')
            while True:
                line = file.readline()
                if '$CHAR_MULT' in line: 
                    break
            line = file.readline()
            words = line.split()
            self.charge = int(words[0])
            self.mult = int(words[1])
        except IOError as err:
            errno, strerror = err.args
            print("I/O error({0}) reading $CHAR_MULT block: {1}".format(errno, strerror))
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise

        # read COORD block
        try:
            while True:
                line = file.readline()
                if '$COORD' in line: 
                    break
            while True:
                line = file.readline()
                if '$END' in line:
                    break
                a = Atom()
                a.read_MKL(line)
                self.atoms.append( a )
        except IOError as err:
            errno, strerror = err.args
            print("I/O error({0}) reading $COORD block: {1}".format(errno, strerror))
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise

        # read BASIS block
        try:
            while True:
                line = file.readline()
                if '$BASIS' in line:
                    break
            line = file.readline()
            for i, a in enumerate(self.atoms):
                if '$END' in line:
                    print("unsufficient entries in $BASIS block")
                    sys.exit(1)
                basis = []
                while True:
                    if '$$' in line or len(line.strip()) == 0:
                        break
                    bf = BasisFunction()
                    line = bf.read_MKL(self.orca_basisorder, line, file)
                    basis.append( bf )
                self.atoms[i].basis_functions = basis
                line = file.readline()
        except IOError as err:
            errno, strerror = err.args
            print("I/O error({0}) reading $BASIS block: {1}".format(errno, strerror))
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise

        # read ALPHA_COORD block
        try:
            while True:
                line = file.readline()
                if '$COEFF_ALPHA' in line:
                    break
            while True:
                line = file.readline()
                if '$END' in line:
                    break
                words = line.split()
                line = file.readline()
                for a in self.atoms:
                    for bf in a.basis_functions:
                        blist = self.orca_basisorder[bf.typ]
                        for bb in blist:                           
                            line = file.readline()
                            clist = [float(w) for w in line.split()]
                            bf.mocoef[bb].extend(clist)
        except IOError as err:
            errno, strerror = err.args
            print("I/O error({0}) reading $COEFF_ALPHA block: {1}".format(errno, strerror))
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise

    def number_of_basis_functions(self):
        n_bf = 0
        for a in self.atoms:
            for bf in a.basis_functions:
                n_bf += len(self.orca_basisorder[bf.typ])
        return n_bf

    def number_of_mos(self):
        n_mo = 0
        d = self.atoms[0].basis_functions[0].mocoef
        l = list(d.keys())
        if len(l) > 0:
            n_mo = len(d[l[0]])    # length of coef list of first key of first basis function
        return n_mo

    def number_of_electrons(self):
        n = 0
        for a in self.atoms:
            n += a.atomic_number
        n -= self.charge
        return n

    def charge_and_mult(self):
        return self.charge, self.mult

    def add_cartesian_coefs(self, transformation_dict):
        n_mo = self.number_of_mos()
        for a in self.atoms:
            for bf in a.basis_functions:
                if bf.typ == 'd' or bf.typ == 'f' or bf.typ == 'g':
                    if 'xx' in bf.mocoef.keys() or 'xxx' in bf.mocoef.keys() or 'xxxx' in bf.mocoef.keys():
                        sys.exit("add_cartesian_coefs: cartesian coefficients exist already")
                    for ss in self.gms_cartorder[bf.typ]:
                        bf.mocoef[ss] = [0.0] * n_mo
                    for bb in self.orca_basisorder[bf.typ]:
                        for cc, ss in zip(transformation_dict[bb], self.gms_cartorder[bf.typ]):
                            for i in range(n_mo):
                                bf.mocoef[ss][i] += cc * bf.mocoef[bb][i]

    def write_amolqc_MOs(self, max_mo=0):
        n_mo = self.number_of_mos()
        if max_mo == 0:
            max_mo = n_mo
        if max_mo > n_mo:
            sys.exit("not enough mos available")

        print(max_mo)
        print(" ")
        for i in range(max_mo):
            bcounter = 0
            print(i + 1)
            for a in self.atoms:
                for bf in a.basis_functions:
                    for ss in self.gms_cartorder[bf.typ]:
                        print("{0:15.8f}".format(bf.mocoef[ss][i]), end='')
                        bcounter += 1
                        if (bcounter%5 == 0):
                            print("")
            if (bcounter%5 != 0):
               print("")



    def write_amolqc_wf(self, basis='gaussian', max_mo=0, title="'generated from moconversion'"):
        print("$general")
        print(" title={0:s}".format(title))
        print(" evfmt=fre, basis={0:s}, jastrow=none".format(basis))
        print(" charge={0:d}, spin={1:d}".format(self.charge, self.mult))
        print("$end")
        print("$geom")
        print(len(self.atoms))
        for a in self.atoms:
            a.write_amolqc()
        print("$end")
        if basis=='gaussian':
            print("$basis")
            for a in self.atoms:
                print(periodic_table[a.atomic_number])
                for bf in a.basis_functions:
                    bf.write_amolqc(sys.stdout)
                print(" ****")
            print("$end")
        print("$mos")
        self.write_amolqc_MOs(max_mo)
        print("$end")
        print("$csfs")
        print("  single restricted")
        print("$end")

if __name__ == "__main__":
    if '-h' in sys.argv or '--help' in sys.argv or len(sys.argv) < 2:
       print('''
    moconvert.py in-file.mkl [max_mo] converts in-file.mkl from orca to amolqc wf form
    writing all or the first max_mo MOs into the wf file.
    ''')
       sys.exit(0)

    fname = sys.argv[1]
    mol = RichOrbitals()
    mol.read_orca_MKL_file(fname)
    mol.add_cartesian_coefs()
    if len(sys.argv) > 2: 
        n_mo = int(sys.argv[2])
        mol.write_amolqc_wf(max_mo=n_mo)
    else:
        mol.write_amolqc_wf()    

