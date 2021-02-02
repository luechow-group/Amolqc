""""
 Copyright (C) 2017-2019 Leonard Reuter

SPDX-License-Identifier: GPL-3.0-or-later
"""

from fractions import Fraction
from .format import format_det_coeff


class Determinant:
    def __init__(self):
        self.coefficient = 1.0  # float
        self.orbital_list = []  # list of integers, beta orbitals as negative ints

    # dets are defined as equal, if their orbital list is equal (independent from coeff)
    def __eq__(self, other):
        if self.__class__.__name__ == other.__class__.__name__:
            if self.orbital_list == other.orbital_list:
                return True
            else:
                return False
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def sort(self, change_sign=True):
        # replacing alpha orbital_list integers by their inverse as a fraction
        for i in range(len(self.orbital_list)):
            if self.orbital_list[i] > 0:
                self.orbital_list[i] = Fraction(1, self.orbital_list[i])
        # bubble sorting
        for i in range(len(self.orbital_list)-1, 0, -1):
            for j in range(i):
                if self.orbital_list[j] < self.orbital_list[j+1]:
                    self.orbital_list[j], self.orbital_list[j+1] = self.orbital_list[j+1], self.orbital_list[j]
                    if change_sign:
                        self.coefficient = -1 * self.coefficient
        # replacing alpha orbital_list fractions with integers again
        for i in range(len(self.orbital_list)):
            if self.orbital_list[i] > 0:
                self.orbital_list[i] = int(Fraction(1, self.orbital_list[i]))

    def write(self, outfile):
        outfile.write(' '+format_det_coeff(self.coefficient))
        for integer in self.orbital_list:
            outfile.write(' '+str(abs(integer)).rjust(2))
        outfile.write('\n')

    def occupation_integer(self):
        occupation_string = ''
        for integer in self.orbital_list:
            occupation_string += str(abs(integer))
        return int(occupation_string)
