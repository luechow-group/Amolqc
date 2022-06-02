""""
 Copyright (C) 2011 Arne Luechow

 SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
if sys.version_info[0] < 3:
    sys.exit('This script requires Python 3')

import numpy as np
from scipy import integrate, linalg
from scipy.special import eval_legendre

"""
This module implements the Gaussian expansion of
STOs using least-squares minimization according to
O-ohata, Taketa, Huzinaga, J. Phys. Soc. Japan 21, 2306 (1966)
1d integration is done numerically (scipy.integrate)
rather than analytically
"""

class Gfctn:
    """Gfctn implements function G_n(alpha) Eq. (2.11)
    using numerical integration"""
    def __init__(self, n, alpha):
        self.alpha = alpha
        self.n = n
    def __call__(self, r):
        return r**self.n * np.exp(-r - self.alpha * r**2)
    def integral(self):
        res, err = integrate.quad(self, 0, np.Inf, limit=100, epsrel=1.e-12, epsabs=1.e-12)
        return res, err

class G1fctn:
    """Gfctn implements function G_n(alpha) Eq. (2.11)
    using numerical integration"""
    def __init__(self, alpha):
        self.alpha = alpha
    def __call__(self, r):
        return (2.0 - r) * np.exp(-0.5*r - self.alpha*r**2) * r**2
    def integral(self):
        res, err = integrate.quad(self, 0, np.Inf, limit=100, epsrel=1.e-12, epsabs=1.e-12)
        return res, err

class AlphaAndCoeffs:
    """AlphaAndCoeffs implements alpha mapping with Legendre polys
       to prevent near linear dependency, G. A. Petersson, JCP 118, 1101, 2003
       constructor: AlphaAndCoeffs(typ,nGTOs) with typ = 1S|2S|2P|3S|3P|3D|4F and
       nGTOs the number of GTOs in expansion"""

    def __init__(self, typ, nGTOs):
        if not nGTOs > 0: raise ValueError("nGTOs must be positive integer")
        self.nGTOs = nGTOs
        if typ == "1S":
            self.nSTO = 1; self.nGTO = 1
        elif typ == "2S":
            self.nSTO = 2; self.nGTO = 1
        elif typ == "2P":
            self.nSTO = 2; self.nGTO = 2
        elif typ == "3S":
            self.nSTO = 3; self.nGTO = 1
        elif typ == "3P":
            self.nSTO = 3; self.nGTO = 2
        elif typ == "3D":
            self.nSTO = 3; self.nGTO = 3
        elif typ == "4F":
            self.nSTO = 4; self.nGTO = 4
        elif typ == "H2s":
            self.nGTO = 1; self.nSTO = 0; self.Hn = 2
        else:
            raise ValueError("AlphaAndCoeffs: unknown value for typ")
        self.alpha = np.zeros(nGTOs, 'd')
        self.c = np.zeros(nGTOs, 'd')
        fak = 1
        for i in range(2,2*self.nSTO+1):
            fak *= i
        self.normSTO = fak**(-0.5) * 2**(self.nSTO + 0.5)
        dfak = 1
        for i in range(3,2*self.nGTO,2):
            dfak *= i 
        # normGTO is only the alpha independent part of N
        self.normGTO = 2.0**(self.nGTO+1) * dfak**(-0.5) * (2.0*np.pi)**(-0.25) 

    def __call__(self, a0):
        for j in range(self.nGTOs):
            lna = 0
            for k in range(len(a0)):
                lna += a0[k] * eval_legendre(k, 2.0 * j / (self.nGTOs - 1.0) - 1.0)
            self.alpha[j] = np.exp(lna)
        V = self.calcCoeffs()
        return V

    def getAlpha(self):
        return self.alpha

    def setAlpha(self, alpha):
        self.alpha = alpha

    def getCoeffs(self):
        return self.c

    def _Pi(self, alp):
        """Pi implements Eq. 2.3 / 2.10 numerically"""
        normGTO = self.normGTO * alp**((2.0*self.nGTO+1.0)/4.0)
        gf = Gfctn(self.nSTO + self.nGTO, alp)
        g,err = gf.integral()
        # if err > 1e-8: raise ValueError
        return self.normSTO * normGTO * g

    def _P1i(self, alp):
        """Pi implements Eq. 2.3 / 2.10 numerically, but with H-like 2s function instead of STO"""
        normGTO = self.normGTO * alp**((2.0*self.nGTO+1.0)/4.0)
        gf = G1fctn(alp)
        g,err = gf.integral()
        # if err > 1e-8: raise ValueError
        return 1.0/np.sqrt(8.0) * normGTO * g

    def calcCoeffs(self):
        """calc Coeffs implements Eq. 2.7/8/9 from
           O-Ohata et al. JPSJ, 21, 2306, 1966"""
        P = np.zeros(len(self.alpha), 'd')
        S = np.zeros((len(self.alpha), len(self.alpha)), 'd')
        if self.nSTO>0:
            for i in range(len(self.alpha)):
                P[i] = self._Pi(self.alpha[i])
                for j in range(i, len(self.alpha)):
                    S[i, j] = ((2 * np.sqrt(self.alpha[i] * self.alpha[j])) / 
                           (self.alpha[i] + self.alpha[j]))**(self.nGTO+0.5)
                    S[j, i] = S[i, j]
        else:
            for i in range(len(self.alpha)):
                P[i] = self._P1i(self.alpha[i])
                for j in range(i, len(self.alpha)):
                    S[i, j] = ((2 * np.sqrt(self.alpha[i] * self.alpha[j])) / 
                           (self.alpha[i] + self.alpha[j]))**(self.nGTO+0.5)
                    S[j, i] = S[i, j]
        self.c = linalg.solve(S, P)
        V = 1 - sum(self.c * P)
        return V

if __name__ == "__main__":
    print(" * * * testing stogto_opt * * *")
    print("testing Gfctn numerical integration")
    g1 = Gfctn(1, 1.0)
    r = 0.0
    print(g1(r))
    print(g1(1.0))
    print(g1.integral())
    
    g1 = Gfctn(2, 2.0)
    print(g1.integral())
    
    print("testing calcCoeffs:")
    # reference: O-Ohata et al. JPSJ, 21, 2306, 1966
    stogto = AlphaAndCoeffs('1S',4)
    alpha = np.array([7.2337, 1.27973, 0.331608, 0.10141])
    stogto.setAlpha(alpha)
    V4 = stogto.calcCoeffs()
    assert abs(V4-6.8e-5) < 1e-6, "V4 1s test failed"
    assert np.linalg.norm(stogto.getAlpha()-alpha) < 1e-6,"getAlpha 1s test failed"
    c = np.array([0.0385,0.2000,0.5204,0.3819])
    assert np.linalg.norm(stogto.getCoeffs()-c) < 1e-3,"calcCoeffs 1s test failed"
    
    
    stogto = AlphaAndCoeffs('1S',10)
    alpha = np.array([1188.35, 156.411, 37.9276, 10.5140, 3.34954, 1.18834, 0.458596,
                   0.191073, 0.0848076, 0.0372356])
    stogto.setAlpha(alpha)
    V10 = stogto.calcCoeffs()
    assert abs(V10-7.6e-8) < 2e-8, "V10 1s test failed"
    c = np.array([8.862e-5,6.228e-4,3.182e-3,1.393e-2,4.961e-2,1.4465e-1,0.3133,0.40100,0.20322,0.015420]) 
    assert np.linalg.norm(stogto.getCoeffs()-c) < 1e-5,"calcCoeffs V10 1s test failed"
    
    # testing 2S expansion into 1S GTOs
    sto2s = AlphaAndCoeffs('2S',4)
    alpha = np.array([11.582,1.9941,0.1611,0.06144])
    c = np.array([-0.012034,-0.05479,0.57810,0.479536])
    sto2s.setAlpha(alpha)
    V4 = sto2s.calcCoeffs()
    assert abs(V4-2.7e-5) < 1e-6,"V4 2s test failed"
    assert np.linalg.norm(sto2s.getCoeffs()-c) < 1e-4,"calcCoeffs V4 2s test failed"
    
    # testing H2s expansion into 1S GTOs
    h2s = AlphaAndCoeffs('H2s',4)
    alpha = np.array([12.0, 2.0, 1.5, 0.1])
    h2s.setAlpha(alpha)
    V4 = h2s.calcCoeffs()
    print(V4,h2s.getCoeffs())
    
    print("all tests passed")


