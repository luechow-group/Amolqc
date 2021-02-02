""""
 Copyright (C) 2017-2019 Leonard Reuter

SPDX-License-Identifier: GPL-3.0-or-later
"""

from .format import format_det_coeff


class Csf:
    def __init__(self):
        self.coefficient = 1.0  # float
        self.determinants = []  # list of determinant objects
        self.occupation = []  # list of strings or ints (different for vb/mo)

    def write(self, outfile):
        outfile.write(format_det_coeff(self.coefficient) + str(len(self.determinants)).rjust(8)+'\n')
        for det in self.determinants:
            det.write(outfile)

    def normalize(self):
        sum_squared_coeffs = 0
        for det in self.determinants:
            sum_squared_coeffs += det.coefficient ** 2
        self.coefficient *= sum_squared_coeffs ** 0.5
        if sum_squared_coeffs != 0.0:
            for det in self.determinants:
                det.coefficient /= sum_squared_coeffs ** 0.5
        else:
            for det in self.determinants:
                det.coefficient = (1.0 / len(self.determinants)) ** 0.5

    def symmetry(self):
        for det in self.determinants:
            pass
