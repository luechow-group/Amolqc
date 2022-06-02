""""
 Copyright (C) 2017-2019 Leonard Reuter

SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
if sys.version_info[0] < 3:
    sys.exit('This script requires Python 3')


class Jastrow:
    def __init__(self):
        self.type = 'none'
        self.lines = []

    def write(self, outfile):
        for line in self.lines:
            outfile.write(line)
