""""
 Copyright (C) 2017-2019 Leonard Reuter

SPDX-License-Identifier: GPL-3.0-or-later
"""


class Jastrow:
    def __init__(self):
        self.type = 'none'
        self.lines = []

    def write(self, outfile):
        for line in self.lines:
            outfile.write(line)
