""""
 Copyright (C) 2017-2019 Leonard Reuter

SPDX-License-Identifier: GPL-3.0-or-later
"""


from .format import format_counter, format_orbital_coefficient


class Orbital:
    def __init__(self):
        self.symmetry = None
        self.coefficients = []  # list of floats

    def write(self, outfile, counter, orbital_format):
        if orbital_format != 'gms':
            outfile.write(str(counter + 1).rjust(4)+'\n')
        for i in range((len(self.coefficients) + 4) // 5):
            if orbital_format == 'gms':
                outfile.write(format_counter(counter + 1)+' '+format_counter(i + 1))
            for j in range(min(5, len(self.coefficients) - 5 * i)):
                outfile.write(format_orbital_coefficient(self.coefficients[i * 5 + j], orbital_format))
            outfile.write('\n')
