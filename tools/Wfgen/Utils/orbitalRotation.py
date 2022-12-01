""""
 Copyright (C) 2018-2019 Leonard Reuter

SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
if sys.version_info[0] < 3:
    sys.exit('This script requires Python 3')


class OrbitalRotation:
    def __init__(self):
        self.orbital_group = []

    def write(self, file):
        for group in self.orbital_group:
            file.write(str(len(group)))
            for index in group:
                file.write(' ' + str(index))
            file.write('\n')
