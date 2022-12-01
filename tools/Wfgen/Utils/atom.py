""""
 Copyright (C) 2017-2019 Leonard Reuter
 Copyright (C) 2020 Arne Luechow

SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
if sys.version_info[0] < 3:
    sys.exit('This script requires Python 3')

from .format import format_coordinate
from .data import periodic_table


class Atom:
    """class containing an atom with coords and possibly basis functions"""   
    def __init__(self):
        self.atomic_number = 0  # integer
        self.position = []
        self.element = ''
        self.charge = 0  # charge for samplin
        self.core_charge = None  # charge different from atomic_number for ECPs
        self.basis_functions = []

    def read_MKL(self, line):
        words = line.split()
        self.atomic_number = int(words[0])
        self.position = [float(words[1]), float(words[2]), float(words[3])]

    def get_element(self):
        self.element = periodic_table[self.atomic_number]  # string
        return self.element

    def write(self, outfile, atomic_charges):
        outfile.write(self.get_element().ljust(2))
        for coordinate in self.position:
            outfile.write(' '+format_coordinate(coordinate))
        if atomic_charges:
            outfile.write(' '+str(self.charge))
        outfile.write('\n')

    def write_molpro(self, outfile, index):
        outfile.write(periodic_table[self.atomic_number].ljust(2)
                      + str(index + 1).rjust(5)
                      + str(self.core_charge).rjust(5)
                      + ' ')
        for coordinate in self.position:
            outfile.write('{: 20.10f}'.format(coordinate))
        outfile.write('\n')
