#
""""
 Copyright (C) 2011 Arne Luechow

 SPDX-License-Identifier: GPL-3.0-or-later
"""
# basisset module: reading and writing basissets
#
from numpy import *
from . import pse
from . import basis

class AtomBasisSet:
    def __init__(self,Z,converter=0,key=0):
        self.Z = Z
        self.basis = []
        self.converter = converter
        self.key = key
    def add(self,bf):
        self.basis.append(bf)
    def read(self,handle):
        try:
            while True:
                bf = basis.readBasisFunction(handle,self.converter,self.key)
                self.add(bf)
        except basis.BasisFunctionError:
            pass
        except:
            raise
    def size(self):
        s = 0
        for bf in self.basis:
            s += bf.nFunctions()
        return s
    def getZ(self):
        return self.Z
    def basisFunction(self,idx):
        return self.basis[idx]
    def changeBasisFunction(self,idx,bf):
        del self.basis[idx]
        self.basis.insert(idx,bf)
    def writeAmolqc(self,handle):
        handle.write(pse.elem(self.Z)+'\n')
        for bf in self.basis:
            bf.writeAmolqc(handle)
        handle.write('****\n')
    def writeGaussian(self,handle):
        handle.write(pse.elem(self.Z)+'  0\n')
        for bf in self.basis:
            bf.writeGaussian(handle)
        handle.write('****\n')
    def writeGamess(self,handle):
        handle.write(pse.elem(self.Z)+'      '+str(self.Z)+'.0    \n')
        for bf in self.basis:
            bf.writeGamess(handle)
        handle.write('****\n')
    def writeGBSGaussian(self,handle):
        handle.write('-'+pse.elem(self.Z)+'  0\n')
        for bf in self.basis:
            bf.writeGaussian(handle)
        handle.write('****\n')

def readGBSFile(filename):
    """read gaussian basis set file (.gbs) and create a list
    mapping element symbols to AtomBasisSet objects"""
    basisset = []
    comments = []
    try: 
        f = open(filename,"r")
        while True:
            line = f.readline()
            if not line or line[0] != '-':
                break
            words = line.split()
            elem = words[0][1:]
            Z = pse.Z(elem)
            atom = AtomBasisSet(Z)
            atom.read(f)
            basisset.append([elem,atom])
        f.close()
    except (IOError,EOFError):
        print("error reading or opening gbs file")
        raise
    return comments,basisset
        
def readEMSLFile(filename):
    """read gaussian basis set file in EMSL format and create a list
    mapping element symbols to AtomBasisset objects"""
    basisset = []
    comments = []
    try: 
        f = open(filename,"r")
        while True:
            line = f.readline()
            if '****' in line:
                break
            comments.append(line)
        while True:
            line = f.readline()
            if not line:
                break
            words = line.split()
            if len(words)<2:
                break
            elem = words[0]
            Z = pse.Z(elem)
            atom = AtomBasisSet(Z)
            atom.read(f)
            basisset.append([elem,atom])
        f.close()
    except (IOError,EOFError):
        print("error reading or opening gbs file")
        raise
    return comments,basisset

def readAmolqcBasisFile(filename,convertFile=None,key=None):
    """read amolqc basis set file  and create a list
    mapping element symbols to AtomBasisset objects
    convertFile contains the data for an optional STOConverter object allowing
    the conversion of STOs into GTOs when required. key is the 
    key used in the convertFile for different expansions"""
    basisset = []
    comments = []
    verbose = False
    if convertFile is not None:
        cvt = basis.STOConverter(convertFile)
    else:
        cvt = None
    try: 
        f = open(filename,"r")
        while True:
            line = f.readline()
            if '****' in line:
                break
            if line[0] == '!':
                comments.append(line)
        while True:
            line = f.readline()
            if not line:
                break
            words = line.split()
            elem = words[0]
            Z = pse.Z(elem)
            atom = AtomBasisSet(Z,cvt,key)
            if verbose: print("read abs file: reading ",elem,Z)
            atom.read(f)
            basisset.append([elem,atom])
        f.close()
    except (IOError,EOFError):
        print("error reading or opening abs file")
        raise
    return comments,basisset

def writeAmolqcBasisFile(basisset,comments,filename):
    """write a basis set in amolqc format"""
    f = open(filename,"w")
    for line in comments:
        f.write(line)
    f.write('****\n')
    for abs in basisset:
        elem,atombf = abs
        atombf.writeAmolqc(f)
    f.close()
    
def writeGBSBasisFile(basisset,filename):
    """write a basis set in gbs format"""
    f = open(filename,"w")
    for abs in basisset:
        elem,atombf = abs
        atombf.writeGBSGaussian(f)
    f.close()
    
def writeGamessBasisFile(basisset,filename):
    """write a basis set in gamess format"""
    f = open(filename,"w")
    for abs in basisset:
        elem,atombf = abs
        atombf.writeGamess(f)
    f.close()
    
def convertBasisFile(intyp,infile,outtyp,outfile,convertFile=None,key=None):
    """convert a basis file in gbs for emsl/g94 format in amolqc format"""
    if intyp == 'gbs':
        comments,basisset = readGBSFile(infile)
    elif intyp == 'emsl':
        comments,basisset = readEMSLFile(infile)
    elif intyp == 'abs':
        comments,basisset = readAmolqcBasisFile(infile,convertFile,key)
    else:
        raise ValueError("convertBasisFile: unknown type of input basis file")
    if outtyp == 'abs':
        writeAmolqcBasisFile(basisset,comments,outfile)
    elif outtyp == 'gbs':
        writeGBSBasisFile(basisset,outfile)
    elif outtyp == 'gms':
        writeGamessBasisFile(basisset,outfile)
    else:
        raise ValueError("convertBasisFile: unknown type of output basis file")
        
    


if __name__ == "__main__":
    print("testing the basisset module")
    f = open("cc-pVTZ-f.gbs")
    while True:
        line = f.readline()
        if '-O' in line:
            break
    OBasis = AtomBasisSet(8)
    OBasis.read(f)
    f.close()
    print('size of basisset:',OBasis.size())
    of = open("test.txt","w")
    OBasis.writeGaussian(of)
    OBasis.writeAmolqc(of)
    bf = basis.STO('1S',8.00)
    OBasis.changeBasisFunction(0,bf)
    OBasis.writeAmolqc(of)
    bf = OBasis.basisFunction(0)
    bf.display()
    bf.setZeta(8.1)
    bf.display()
    OBasis.writeAmolqc(of)
    
    f = open("cc-pVTZ-f.gbs")
    while True:
        line = f.readline()
        if '-H' in line:
            break
    HBasis = AtomBasisSet(1)
    HBasis.read(f)
    f.close()
    print('size of basis set:',HBasis.size())
    HBasis.writeGaussian(of)
    HBasis.writeAmolqc(of)
    print("type of AtomBasisSet:",type(HBasis))
    comments,b6311 = readGBSFile("6311.gbs")
    print("size of basis set:",len(b6311))
    cbf = b6311[5]
    print("6-th element ",cbf[0],pse.Z(cbf[0]))
    cbf[1].writeAmolqc(of)
    of.close()
    print("converting 6-311Gdp to amolqc abs format")
    convertBasisFile('emsl','6311gdp.g94','abs','6-311Gdp.abs')
    print("reading amolqc abs file")
    comments,s311 = readAmolqcBasisFile('6-311Gdp.abs')
    # init part
    bfList = [ ['H',0],['C',0] ]   # list containing element and index of basis function for element
    keys = []
    for bf in s311:
        keys.append(bf[0])
    ptrList = []
    for bf in bfList:
        elem,k = bf
        idx = keys.index(elem)
        atombs = s311[idx][1]
        ptrList.append([atombs,k])
    # function call part
    a = [1.1,6.1]
    for p,z in zip(ptrList,a):
        atombs,idx = p
        sto = basis.STO('1S',z)
        atombs.changeBasisFunction(idx,sto)
    writeAmolqcBasisFile(s311,[],'test.abs')       
    
    
    
    
    
    
