""""
 Copyright (C) 2020 Arne Luechow

SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
if sys.version_info[0] < 3:
    sys.exit('This script requires Python 3')

class BasisFunction:
    """class containing a basis function (with p, d degeneracy) and MO coefficients"""

    def __init__(self):
        self.typ = ''
        self.coef = 1.0
        self.cntrctn = []
        self.mocoef = {}

    def read_MKL(self, basisorder, line, file):
        words = line.split()
        self.typ = words[1].lower()
        self.coef = float(words[2])
        self.cntrct = []
        self.mocoef = {}
        line = file.readline()
        while True:
            if line[0:2] != '  ':
                break
            words = line.split()
            self.cntrct.append( [ float(words[0]), float(words[1]) ] ) 
            for bb in basisorder[self.typ]:
                self.mocoef[bb] = [] 
            line = file.readline()
        return line

    def write_amolqc(self, out):
        out.write("{0:s}   {1:3d} {2:12.8f}\n".format(str.upper(self.typ), len(self.cntrct),  self.coef))
        for cntrct in self.cntrct:
            out.write("  {0:12.5f}  {1:12.8f}\n".format(cntrct[0], cntrct[1]))

