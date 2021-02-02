#!/usr/bin/env python
""""
 Copyright (C) 2011 Arne Luechow

 SPDX-License-Identifier: GPL-3.0-or-later
"""
# basis module: convert formats and STO/GTO
#
from numpy import *
import re


class BasisFunctionError(Exception):
    def __init__(self,value="no basis function"):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Basis:
    """ base class for basis functions"""
    def __init__(self,typ):
        self.typ = typ
    def display(self):
        print("Basis:",str(self.typ))
    def nFunctions(self):
        return 1
    def basisType(self):
        return self.typ
    def writeAmolqc(self,handle): pass
    def writeGaussian(self,handle): pass

class CGTO(Basis):
    """CGTO: contracted GTOs. constructor and write methods for
    output, amolqc and gaussian formats. This type allows for Pople
    SP functions that share alphas for s and p functions
    CGTO type are S,P,SP,D,F"""

    def __init__(self,typ,alpha,coeffs,coeffp=[],ccc=[]):
        """ arguments: string typ, list alpha, list coeffs, optional list coeffp
        optional cusp correction coefficients ccc"""
        assert len(alpha)>0 and len(alpha)==len(coeffs),"CGTO: equal sizes required"
        assert typ in ['S','P','SP','D','F','G'],"CGTO: unknown type"
        if typ=='SP':
            assert coeffp != [],"CGTO: SP functions require 2nd coeffs (for p)"
        Basis.__init__(self,typ)
        self.alpha = list(alpha)
        self.c = list(coeffs)
        self.cp = list(coeffp)
        self.ccc = list(ccc)
        self.nGTO = len(alpha)

    def display(self):
        print("CGTO:",self.typ, self.nGTO)
        if self.typ != 'SP':
            for a,c in zip(self.alpha,self.c):
                print(a,c)
        else:
            for a,cs,cp in zip(self.alpha,self.c,self.cp):
                print(a,cs,cp)

    def getnGTO(self):
        """ number of primitive GTOs"""
        return self.nGTO

    def nFunctions(self):
        """ number of basis functions (2 for Pople SP type)"""
        if self.typ == 'SP':
            n = 2
        else:
            n = 1
        return n

    def writeAmolqc(self,handle):
        """ write basis function in amolqc format"""
        if self.typ != 'SP':
            ccstr = ''
            if self.ccc != []:
                for a in self.ccc:
                    ccstr += " "+str(a)
            handle.write(self.typ+" GTO "+str(self.nGTO)+"    "+ccstr+"\n")
            for alp,c in zip(self.alpha,self.c):
                s = ' {0:20.12e} {1:20.12e}'.format(alp,c)
                s1 = s.replace('e','D')
                handle.write(s1+"\n")
        else:
            handle.write("S GTO "+str(self.nGTO)+"\n")
            for alp,c in zip(self.alpha,self.c):
                s = ' {0:20.12e} {1:20.12e}'.format(alp,c)
                s1 = s.replace('e','D')
                handle.write(s1+"\n")
            handle.write("P GTO "+str(self.nGTO)+"\n")
            for alp,c in zip(self.alpha,self.cp):
                s = ' {0:20.12e} {1:20.12e}'.format(alp,c)
                s1 = s.replace('e','D')
                handle.write(s1+"\n")

    def writeGaussian(self,handle):
        ccstr = ''
        if self.ccc != []:
            for a in self.ccc:
                ccstr += " "+str(a)
        handle.write(" "+self.typ+"  "+str(self.nGTO)+" 1.00  "+ccstr+"\n")
        if self.typ != 'SP':
            for alp,c in zip(self.alpha,self.c):
                s = ' {0:20.12e} {1:20.12e}'.format(alp,c)
                s1 = s.replace('e','D')
                handle.write(s1+"\n")
        else:
            for alp,c,cp in zip(self.alpha,self.c,self.cp):
                s = ' {0:20.12e} {1:20.12e} {2:20.12e}'.format(alp,c,cp)
                s1 = s.replace('e','D')
                handle.write(s1+"\n")

    def writeGamess(self,handle):
        ccstr = ''
        if self.ccc != []:
            for a in self.ccc:
                ccstr += " "+str(a)
        handle.write(" "+str(self.typ)+"     "+str(self.nGTO)+"\n")
        if self.typ != 'SP':
            i = 0
            for alp,c in zip(self.alpha,self.c):
                i += 1
                s = ' {0:6d} {1:20.10f} {2:14.10f}'.format(i,alp,c)
                handle.write(s+"\n")
        else:
            print("writeGamess: SP basis functions not implemented")
            raise


class STOConverter:
    """ STOConverter: converts an STO to GTOs using expansion data from file"""
    def __init__(self,filename):
        self.d = { '1S':{}, '2S':{}, '2P':{}, '3S':{}, '3P':{}, '3D':{} , '4F':{} }
        try:
            f = open(filename,"r")
        except IOError as e:
            print("STOConverter: reading convert file failed with: ",e)
            raise
        reg = re.compile("1S|2S|2P|3S|3P|3D|4F")
        while True:
            line = f.readline()
            if "****" in line:
                break
            m = reg.match(line)
            if m:
                words = line.split()
                dict = self.d[words[0]]
                key = words[1]
                words = f.readline().split()
                nGTO = int(words[0])
                alpha = zeros(nGTO,'d')
                coeffs = zeros(nGTO,'d')
                for k in range(nGTO):
                    words = f.readline().split()
                    alpha[k] = float(words[0])
                    coeffs[k] = float(words[1])
                typ = m.group(0)[1:]
                cgto = CGTO(typ,alpha,coeffs)
                dict[key] = cgto
        f.close()

    def display(self):
        for typ in self.d.keys():
            dict = self.d[typ]
            for key in dict.keys():
                print(key,dict[key].getnGTO())

    def toCGTO(self,typ,zeta,key=None):
        if typ not in ["1S","2S","2P","3S","3P","3D","4F"]:
            raise ValueError("toCGTO: unknown typ")
        dict = self.d[typ]
        keys = [item[0] for item in dict.items()]
        if key is None:
            if len(keys)==0: raise ValueError("toCGTO: no expansion found")
            key = keys[0]
        else:
            if key not in keys: raise ValueError("toCGTO: key not found")
        cgto = dict[key]
        # now use scaling theorem to get STO expansion for zeta
        # use array not list
        alpha = array(cgto.alpha)
        alpha *= zeta**2
        cgto1 = CGTO(typ[1:],alpha,cgto.c)
        return cgto1

class STO(Basis):
    """STO with 'typ' and 'zeta'. 'typ' is only 1S,2P,3D,4F
    converter and key handles allow conversion to CGTO"""

    def __init__(self,typ,zeta,converter=None,key=None):
        if typ not in ['1S','2S','2P','3S','3P','3D','4F']:
            raise ValueError("STO: unknown type")
        Basis.__init__(self,typ)
        self.zeta = zeta
        self.converter = converter
        self.key = key
        self.cgto = 0

    def display(self):
        print("STO:",self.typ, self.zeta)

    def setZeta(self,zeta):
        self.zeta = zeta

    def writeAmolqc(self,handle):
        handle.write(self.typ+" STO "+str(self.zeta)+"\n")

    def writeGaussian(self,handle):
        if self.converter != 0:
            if self.key == 0:
                self.cgto = self.converter.toCGTO(self.typ,self.zeta)
            else:
                self.cgto = self.converter.toCGTO(self.typ,self.zeta,self.key)
            self.cgto.writeGaussian(handle)
        else:
            raise BasisFunctionError("Gaussian output requires STO-GTO converter")

def readBasisFunction(handle,converter=None,key=None):
    """readBasisFunction: utility function to read basis function from file
    both amolqc and gaussian formats are recognized. converter/key can be supplied
    to STO functions for later conversion to GTOs. The basis function is returned"""

    line = handle.readline()
    if '****' in line:
        raise BasisFunctionError("not a basis function")
    elif 'STO' in line:
        words = line.split()
        typ = words[0]
        zeta = float(words[2])
        bf = STO(typ,zeta,converter,key)
    elif 'GTO' in line:
        words = line.split()
        typ = words[0]
        nGTO = int(words[2])
        if len(words) > 3:
            assert len(words)==7,"readBasisFunction: cusp correction incomplete"
            ccc = [float(words[3]),float(words[4]),float(words[5]),float(words[6])]
        else:
            ccc = []
        alpha = [0.0]*nGTO
        coeffs= [0.0]*nGTO
        for k in range(nGTO):
            words = handle.readline().split()
            alpha[k] = float(words[0].replace('D','e'))
            coeffs[k] = float(words[1].replace('D','e'))
        bf = CGTO(typ,alpha,coeffs,ccc=ccc)
    else:
        words = line.split()
        typ = words[0]
        nGTO = int(words[1])
        if len(words) > 3:
            assert len(words)==7,"readBasisFunction: cusp correction incomplete"
            ccc = [float(words[3]),float(words[4]),float(words[5]),float(words[6])]
        else:
            ccc = []
        alpha = [0.0]*nGTO
        coeffs = [0.0]*nGTO
        if typ == 'SP':
            coeffp = [0.0]*nGTO
            for k in range(nGTO):
                words = handle.readline().split()
                alpha[k] = float(words[0].replace('D','e'))
                coeffs[k] = float(words[1].replace('D','e'))
                coeffp[k] = float(words[2].replace('D','e'))
            bf = CGTO(typ,alpha,coeffs,coeffp,ccc=ccc)
        else:
            for k in range(nGTO):
                words = handle.readline().split()
                alpha[k] = float(words[0].replace('D','e'))
                coeffs[k] = float(words[1].replace('D','e'))
            bf = CGTO(typ,alpha,coeffs,ccc=ccc)
    return bf


if __name__ == "__main__":
    print("testing the basis module")
    b = Basis("1S")
    b.display()
    alpha = [100.0, 10.0, 1.0, 0.1]
    coeffs = [0.01, 0.3, 0.5, 0.1]
    nGTO = len(alpha)
    g1 = CGTO("S",alpha,coeffs)
    g1.display()
    assert g1.getnGTO()==4,"getnGTO failed"
    zeta = 6.5
    s1 = STO("1S",zeta)
    s1.display()
    cvt = STOConverter("sto2gto.dat")
    cvt.display()
    s2 = STO("1S",zeta,converter=cvt,key="OPT10")
    s2.display
    f = open("test.out","w")
    g1.writeAmolqc(f)
    g1.writeGaussian(f)
    s1.writeAmolqc(f)
    s2.writeAmolqc(f)
    s2.writeGaussian(f)
    zeta = 7.0
    s2.setZeta(zeta)
    s2.writeAmolqc(f)
    #s2.writeGaussian(f)
    fi = open("cc-pVTZ-f.gbs","r")
    line = fi.readline()
    bf = readBasisFunction(fi)
    bf.display()
    f.write("gto from cc-pCTZ-f\n")
    bf.writeAmolqc(f)
    f.close()




