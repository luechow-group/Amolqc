#
""""
 Copyright (C) 2011 Arne Luechow

 SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
if sys.version_info[0] < 3:
    sys.exit('This script requires Python 3')

import numpy as np
from scipy.optimize import optimize
from Basis import stogto_opt

def getGTOexpansion(alpha,typ):
    """getGTOexpansion(alpha,typ): calculates coefficients and least-squares variance
    for given alphas-array for STO typ='1S'|'2S'|'2P'|'3S'|'3P'|'3D' """
    assert typ in ['1S','2S','2P','3S','3P','3D'], "runSTOGTOOptimization: illegal argument"
    #assert type(alphas)==type(np.array[0.0]), "runSTOGTOOptimization: illegal array argument"
    nGTOs = len(alpha)
    sto = stogto_opt.AlphaAndCoeffs(typ,nGTOs)
    sto.setAlpha(alpha)
    V = sto.calcCoeffs()
    c = sto.getCoeffs()
    print("V = ",V)
    print("printing i alpha_i coeff_i")
    print(typ+" INP"+str(nGTOs)+" (zeta=1.0) expansion from stogto_opt.py")
    print(len(alpha))
    for i in range(len(alpha)):
        print(' {0:20.10f} {1:20.10f}'.format(alpha[i], c[i]))
    print(" ")
    l = []
    for i,j in zip(alpha,c):
        l.append([i,j])
    return l

def runSTOGTOOptimization(a0,typ,nGTOs):
    """runSTOGTOOptimization(a0,typ,nGTOs)
    optimizes the expansion of an STO of typ='1S'|'2S'|'2P'|'3S'|'3P'|'3D' into nGTOs
    the vector a0 contains e.g. 4 floats that expand into nGTOs alphas avoiding
    near linar dependency (algorithm of Petersson)"""
     
    assert typ in ['1S','2S','2P','3S','3P','3D','4F','H2s'], "runSTOGTOOptimization: illegal argument"
    #assert type(a0)==type(np.array[0.0]), "runSTOGTOOptimization: illegal array argument"
    assert nGTOs > 0, "runSTOGTOOptimization: nGTOs required"
    stoopt = stogto_opt.AlphaAndCoeffs(typ,nGTOs)
    print(" start values a0=", a0)
    V = stoopt(a0)
    print("start variance V0=", V)
    print("optimizing expansion into ",nGTOs," GTOs...")
    a0opt = optimize.fmin(stoopt, a0, xtol=1.e-14)
    print("done")
    print(" final A_k coeffs: ", a0opt)
    print(" final lsq value V"+str(nGTOs)+" = ", stoopt(a0opt))
    a = stoopt.getAlpha()
    c = stoopt.getCoeffs()
    print(typ+" OPT"+str(nGTOs)+" (zeta=1.0) expansion from stogto_opt.py")
    print(len(a))
    for i in range(len(a)):
        print(' {0:20.10e} {1:20.10e}'.format(a[i], c[i]))
    print(" ")
    l = []
    for i,j in zip(a,c):
        l.append([i,j])
    return l

# reproduce O-Ohata et al. results
alpha = np.array([11.582,1.9941,0.1611,0.06144])
getGTOexpansion(alpha,'2S')

alpha = np.array([80.1823, 15.1962, 4.17165, 1.36064,0.207424, 0.0973524, 0.0558673, 0.0382022])
getGTOexpansion(alpha,'2S')

alpha = np.array([432.099, 60.6176, 14.1623, 3.98351, 1.3073, 0.463514, 0.201066, 0.0991734, 0.049785, 0.0214081])
getGTOexpansion(alpha,'2S')

alpha = np.array([1.46566, 0.384173, 0.139164, 0.0577971])
getGTOexpansion(alpha,'2P')

alpha = np.array([39.9819, 11.5715, 3.94709, 1.42281, 0.548285, 0.227070, 0.100995, 0.0467295])
getGTOexpansion(alpha,'2P')

alpha = np.array([1.85790, 0.191764, 0.0872582, 0.0424396])
getGTOexpansion(alpha,'3P')

alpha = np.array([1.28394, 0.387764, 0.147129, 0.061157])
getGTOexpansion(alpha,'3D')


a0 = np.array([-1.0, -3.5, 0.1, -0.05])
sto = stogto_opt.AlphaAndCoeffs('1S',10)
print(sto(a0))
print(sto.getAlpha())

#a0 = np.array([-1.0, -3.5, 0.5, -0.1])
#l = runSTOGTOOptimization(a0,'H2s',14)

print(" ************ starting new optimization: **********")
l = runSTOGTOOptimization(a0,'4F',14)
###print(l

    