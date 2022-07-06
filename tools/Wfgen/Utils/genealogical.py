#!/usr/bin/env python

"""
 Copyright (C) 2017 Arne Luechow

SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
if sys.version_info[0] < 3:
   sys.exit('This script requires Python 3')

from fractions import Fraction
from math import sqrt

def P(i, j, b):
   # permute bits i and j in binary number b
   bnew = b
   # read bit i
   imask = 1 << i
   bi = ( b & imask != 0 )
   jmask = 1 << j
   bj = ( b & jmask != 0 )
   if bi:
      bnew = bnew | jmask
   else:
      bnew = bnew & ~ jmask
   if bj:
      bnew = bnew | imask
   else:
      bnew = bnew & ~ imask
   return bnew

class spinfunction:
   # spin representation as bit alpha==0, beta==1
   # primitive spin function as binary number
   # N: number of spins (electrons)
   # general spin is represented as coefficient vector (c)
   # i.e. linear combination of primitive spin functions
   # i.e. element of vector space V_N (size 2**N) spanned by primitive spin functions
   # the vector index represents the primitive spin function
   # optional: gen representing the genealogy as string
   def __init__(self, N, gen = None, b = None):
      self.N = N
      self.c = [0 for i in range(2**N)]
      if b is not None:
         self.c[b] = 1
      self.gen = ""
      if gen is not None:
         self.gen = gen

   def __str__(self):
      s = self.gen
      for k in range(2**self.N):
         if self.c[k] != 0:
            ss = ""
            mask = 1
            for i in range(self.N):
               if (mask & k) == 0:
                  ss = "a" + ss
               else:
                  ss = "b" + ss
               mask <<= 1   
            s += " " + str(self.c[k]) + "*" + ss
      return s

   def __add__(self, other):
      assert(self.N == other.N)
      sfnew = spinfunction(self.N, self.gen)
      for k in range(2**self.N):
         sfnew.c[k] = self.c[k] + other.c[k]
      return sfnew

   def __rmul__(self, scal):
      sfnew = spinfunction(self.N, self.gen)
      for k in range(2**self.N):
         sfnew.c[k] = scal * self.c[k]
      return sfnew

   def append_a(self):
      sfnew = spinfunction(self.N + 1, self.gen)
      for k in range(2**self.N):
         sfnew.c[2*k] = self.c[k]
      return sfnew

   def append_b(self):
      sfnew = spinfunction(self.N + 1, self.gen)
      for k in range(2**self.N):
         sfnew.c[2*k + 1] = self.c[k]
      return sfnew

   def szi(self, i):
      # s_z(i) operator
      assert(i < self.N)
      sfnew = spinfunction(self.N, self.gen)
      b = 2**i
      for k in range(2**self.N):
         if self.c[k] != 0:
            # test i-th bit
            if k & b != 0:
               sfnew.c[k] = - Fraction(1, 2) * self.c[k]
            else:
               sfnew.c[k] =   Fraction(1, 2) * self.c[k]
      return sfnew

   def s2i(self, i):
      # s^2(i) operator
      assert(i < self.N)
      sfnew = spinfunction(self.N, self.gen)
      for k in range(2**self.N):
         if self.c[k] != 0:
            sfnew.c[k] = Fraction(3, 4) * self.c[k]
      return sfnew 

   def splusi(self, i):
      assert(i < self.N)
      sfnew = spinfunction(self.N, self.gen)
      b = 1 << i
      mask = ~ b      # complement
      for k in range(2**self.N):
         if self.c[k] != 0:
            if k & b != 0:
               knew = k & mask
               sfnew.c[knew] += self.c[k]
      return sfnew

   def sminusi(self, i):
      assert(i < self.N)
      sfnew = spinfunction(self.N, self.gen)
      b = 2**i
      for k in range(2**self.N):
         if self.c[k] != 0:
            if k & b == 0:
               knew = k | b
               sfnew.c[knew] += self.c[k]
      return sfnew

   def s2(self):
      # implementing Dirac's formula
      # only for primitive spin functions, here unit vector!
      count = 0 
      for k in range(2**self.N):
         if self.c[k] != 0:
            if self.c[k] == 1:
               count = count + 1
               k0 = k
            else:
               count = count + 2
      assert(count == 1)
      sfnew = spinfunction(self.N, self.gen)
      sfnew.c[k0] = - Fraction( self.N * (self.N - 4), 4) 
      for i in range(self.N):
         for j in range(i+1, self.N):
            kk = P(i, j, k0)
            sfnew.c[kk] = 1
      return sfnew

def inner_prod(sf1, sf2):
   # note that coefficients have correct sign but squared absolute value
   assert(sf1.N == sf2.N)
   res = 0
   for k in range(2**sf1.N):
      if sf1.c[k] != 0 and sf2.c[k] != 0:
         prod = sf1.c[k] * sf2.c[k]
         if (prod > 0):
            res = res + sqrt(prod)
         else:
            res = res - sqrt(-prod)
   return res

def create_genealogical_spin_functions(Nmax, singlet = None):
   # S1: S if integer, S - 1/2 if S half-integer
   # M1: M + S (always integer >= 0!)
   # first set up X(N, S1, M) as empty k list
   # optional: singlet=True: calculate only those spin function required for 
   # X(Nmax, S=0, M=0)
   X = []
   for N in range(Nmax + 1):
      Slist = []
      for S1 in range(N//2 + 1):      # new integer division in Python 3
         if N%2 == 0:
            SMax = 2*S1 + 1
         else:
            SMax = 2*S1 + 2
         Mlist = [[] for M1 in range(SMax)]   # empty k list
         Slist.append(Mlist)
      X.append(Slist)

   # append |a> and |b> to empty k lists as X(1, 0, 1, 0) and X(1, 0, 0, 0)
   N = 1
   sa = spinfunction(N, gen="/", b=0)   # |a>
   sb = spinfunction(N, gen="/", b=1)   # |b>
   X[1][0][0].append(sb)
   X[1][0][1].append(sa)

   # now genealogical construction (following Pauncz book)
   for N in range(2, Nmax + 1):
      if singlet is None:
         S0Max = N//2          # new integer division in Python 3
      elif singlet:
         S0Max = min(N//2, (Nmax - N)//2)
      for S1 in range(S0Max + 1):      
         if N%2 == 0:
            SMax = 2*S1 + 1
            S = S1
            Splushalf = S1
            Sminushalf = S1 - 1
         else:
            SMax = 2*S1 + 2
            S = S1 + Fraction(1, 2)
            Splushalf = S1 + 1
            Sminushalf = S1
         for M1 in range(SMax):
            M = M1 - S
            # subtraction if possible
            if S1 < N//2:
               Mminushalf = int(M - Fraction(1, 2) + S + Fraction(1, 2))
               Mplushalf = int(M + Fraction(1, 2) + S + Fraction(1, 2))
               sf1list = []
               sf2list = []
               if Mminushalf >= 0:
                  for k in range(len(X[N - 1][Splushalf][Mminushalf])):
                     sfnew = Fraction(-(S - M + 1), 2*S + 2) * X[N - 1][Splushalf][Mminushalf][k].append_a()
                     sf1list.append(sfnew)
               if Mplushalf < len(X[N - 1][Splushalf]):
                  for k in range(len(X[N - 1][Splushalf][Mplushalf])):
                     sfnew = Fraction(S + M + 1, 2*S + 2) * X[N - 1][Splushalf][Mplushalf][k].append_b()
                     sf2list.append(sfnew)
               #print(N, S, M, "sub:")
               if len(sf1list) == 0:
                  for sf2 in sf2list:
                     sf2.gen = sf2.gen + "\\"
                     #print(sf2)
                     X[N][S1][M1].append(sf2)
               elif len(sf2list) == 0:
                  for sf1 in sf1list:
                     sf1.gen = sf1.gen + "\\"
                     #print(sf1)
                     X[N][S1][M1].append(sf1)
               else:
                  for sf1 in sf1list:
                     for sf2 in sf2list:
                        if sf1.gen == sf2.gen:
                           sfnew = sf1 + sf2
                           sfnew.gen = sf1.gen + "\\"
                           #print(sfnew)
                           X[N][S1][M1].append(sfnew)
            # addition if possible (if S (not S1) > 0)
            if S > 0:
               Mminushalf = int(M - Fraction(1, 2) + S - Fraction(1, 2))
               Mplushalf = int(M + Fraction(1, 2) + S - Fraction(1, 2))
               sf1list = []
               sf2list = []               
               if Mminushalf >= 0:
                  for k in range(len(X[N - 1][Sminushalf][Mminushalf])):
                     sfnew = Fraction(S + M, 2*S) * X[N - 1][Sminushalf][Mminushalf][k].append_a()
                     sf1list.append(sfnew)
               if Mplushalf < len(X[N - 1][Sminushalf]):
                  for k in range(len(X[N - 1][Sminushalf][Mplushalf])):
                     sfnew = Fraction(S - M, 2*S) * X[N - 1][Sminushalf][Mplushalf][k].append_b()
                     sf2list.append(sfnew)
               #print(N, S, M, "add:")
               if len(sf1list) == 0:
                  for sf2 in sf2list:
                     sf2.gen = sf2.gen + "/"
                     #print(sf2)
                     X[N][S1][M1].append(sf2)
               elif len(sf2list) == 0:
                  for sf1 in sf1list:
                     sf1.gen = sf1.gen + "/"
                     #print(sf1)
                     X[N][S1][M1].append(sf1)
               else:
                  for sf1 in sf1list:
                     for sf2 in sf2list:
                        if sf1.gen == sf2.gen:
                           sfnew = sf1 + sf2
                           sfnew.gen = sf1.gen + "/"
                           #print(sfnew)
                           X[N][S1][M1].append(sfnew)
   return X

def printX(X, N, S, M):
   # see definition of S1, M1 in class spinfunction
   if S == int(S):
      S1 = int(S)
   else:
      S1 = int(S - Fraction(1, 2))
   M1 = int(S + M)
   Xlist = X[N][S1][M1]
   for sf in Xlist:
      print(sf)
