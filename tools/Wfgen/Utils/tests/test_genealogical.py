"""
 Copyright (C) 2017 Arne Luechow
 Copyright (C) 2019 Leonard Reuter

SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
if sys.version_info[0] < 3:
    sys.exit('This script requires Python 3')

from Wfgen.Utils.genealogical import *
import unittest as ut


class TestGenealogical(ut.TestCase):
    def test_create_genealogical_spin_functions(self):
        X = create_genealogical_spin_functions(4, singlet=True)
        self.assertEqual(len(X[4][0][0]),2)
        self.assertEqual(inner_prod(X[4][0][0][0],X[4][0][0][1]),0.0)
        self.assertEqual(inner_prod(X[4][0][0][0],X[4][0][0][0]),1.0)
        self.assertEqual(inner_prod(X[4][0][0][1],X[4][0][0][1]),1.0)

if __name__ == '__main__':
    ut.main()
