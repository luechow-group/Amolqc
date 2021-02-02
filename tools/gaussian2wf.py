#!/usr/bin/env python

""""
 TODO UNKNOWN ORIGIN AND AUTHOR (OLDER THAN 2011?)
 Copyright (C) 2011-2012 Arne Luechow
 Copyright (C) 2013 Alexander Sturm
 Copyright (C) 2015-2016 Christoph Schulte
 Copyright (C) 2018 Leonard Reuter

 SPDX-License-Identifier: GPL-3.0-or-later
"""
import sys,getopt
#
def get1st(item):
  return item[0]
#
def gaussian2wf(fname,basis,method='R',fortfilename='fort.7',title='no title',
    extension='.log',stdorient=True,coeffcutoff=0,dets=True,ecp='none',ecpatoms=[]):
    """gaussian2wf constructs initial amolqc wf file using 'fname'.log, the 'basis' name
    The orbitals are read from the Gaussian punch file 'fortfilename'
    """
    if not(ecp =='none'):
        aebasis=basis
        basis='diff'

    logfilename = fname+extension
    wffilename = fname+'.wf'
    if stdorient:
        orient = 'Standard orientation'
    else:
        orient = 'Input orientation'

    if ("casscf" in method) or ("CASSCF" in method):
        active_space=method[method.find('(')+1:method.find(')')]
        orbitals=active_space.split(',')
        casorb = int(orbitals[1])
        caselec = int(orbitals[0])
        if method[0] != 'R' and method[0] != 'U':
            method='R-'+method
    else:
        if not (method[0]=='R' or method[0]=='U'):
            print("gaussian2wf requires method starting with 'R' or 'U'")
            sys.exit(1)

    print(' gaussian2wf:')
    print(' reading ',logfilename,'  and ',fortfilename)
    print(' writing ',wffilename)
    print(' with basis ',basis,'  and title ',title)
    print(' method ',method,'  and extension ',extension)
    print(' orient ',orient)

    try:
        logfile = open(logfilename,'r')
    except IOError as e:
        print('Could not open logfile ',logfilename)
        print(e)
        raise
    try:
        fortfile = open(fortfilename,'r')
    except IOError as e:
        print('Could not open punch file ',fortfilename)
        print(e)
        raise

    wffile = open(wffilename,'w')
    pse = ['X','H','He',
       'Li','Be','B','C','N','O','F','Ne',
       'Na','Mg','Al','Si','P','S','Cl','Ar',
       'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr']
    electronegativity = [float('nan'),2.2,3.5,
        0.98,1.57,2.04,2.55,3.04,3.44,3.98,4.3,
        0.93,1.31,1.61,1.9,2.19,2.58,3.16,3.6,
        0.82,1.0,1.36,1.54,1.63,1.66,1.55,1.83,1.88,1.91,1.9,1.65,1.81,2.01,2.18,2.55,2.96,3.2]

    while True:
        line = logfile.readline()
        if not line: raise EOFError("Charge not found in log file")
        if 'Charge =' in line:
            words = line.split()
            charge = int(words[2])
            mult = words[5]
            break
    print(' charge=',str(charge), 'mult=',mult)

    geometry = readGeometry(logfile, orient, pse)

    while True:
        line = logfile.readline()
        if not line: raise EOFError('number of alpha/beta electrons not found in log file')
        if 'alpha electrons' in line:
            words = line.split()
            nalpha = int(words[0])
            nbeta = int(words[3])
            break
    logfile.close()

    #reopen log-file to find dets
    findDets = dets and (("casscf" in method) or ("CASSCF" in method))
    if findDets:
        try:
            logfile = open(logfilename,'r')
        except IOError as e:
            print('Could not open logfile ',logfile)
            print(e)
            raise
        configs_read=[]
        eigenvectors=[]
        read=False
        while True:
            line = logfile.readline()
            if not line: raise EOFError('CAS - configuration not found')
            if 'Configuration' in line:
                while ('Symmetry' in line) and ('Configuration' in line):
                    words=line.split()
                    configs_read.append(words[4])
                    line = logfile.readline()
            if 'EIGENVALUES AND  EIGENVECTORS' in line:
                 eigenvectors=[0.0]*len(configs_read)
                 for i in range(5):
                    line = logfile.readline()
                 while not 'Final' in line:
                     words=line.split('(')
                     words.sort()  # sort list
                     for word in words:
                         if ')' in word:
                            temp=word.split(')')
                            eigenvectors[int(temp[0])-1]=(float(temp[1]))
                     line=logfile.readline()
                 read=True
            if read: break
        #print(configs_read)
        #print(eigenvectors)
        logfile.close()
     # end if  '  casscf in method '

    print("log-File read")

    # number of orbitals
    nalphaorb = nalpha
    nbetaorb = nbeta

    #building dets
    if findDets:
        det=[]
        excess_alpha = int(mult) - 1
        core_alpha = nalpha - (caselec + excess_alpha) // 2
        core_beta = nbeta - (caselec - excess_alpha) // 2
        nalphaorb = core_alpha + casorb
        nbetaorb = core_beta + casorb

        for config in configs_read:
           det_alpha=''
           det_beta=''
           for i in range(core_alpha):
               det_alpha = det_alpha +' '+ str(i+1)
           for i in range(core_beta):
               det_beta = det_beta +' '+ str(i+1)

           for i in range(len(config)):
               if config[i] == '1':
                   det_alpha = det_alpha +' '+ str(i+1 + core_alpha)
                   det_beta = det_beta +' '+ str(i+1 + core_beta)
               elif config[i]=='a':
                  det_alpha = det_alpha +' '+ str(i+1 + core_alpha)
               elif config[i]=='b':
                  det_beta = det_beta +' '+ str(i+1 + core_beta)
           det.append(det_alpha + det_beta)

    # reading fort.7
    # alpha block
    na = 0
    alphamos = []
    line = fortfile.readline()

    while line and 'Alpha' not in line:
        line = fortfile.readline()
    if not line:
        sys.exit("ERROR: wrong file format of fort.7")

    while na < nalphaorb:
        if 'Alpha' in line:
            words = line.split()
            idx = int(words[0])
            mosblock = []
            mosblock.append(line)
            line = fortfile.readline()
            while 'MO' not in line:
                mosblock.append(line)
                line = fortfile.readline()
            na += 1
            alphamos.append(mosblock)
        elif 'Beta' in line:
            sys.exit("ERROR: not enough alpha MOs in fort.7")
        else:
            sys.exit("ERROR: wrong file format of fort.7")

    if method[0] == 'U':
        nb = 0
        betamos = []
        if nbeta > 0:
            while 'Beta' not in line:
                line = fortfile.readline()
                if not line: raise EOFError("no beta orbitals found in punch file")
        while nb < nbetaorb:
            words = line.split()
            idx = int(words[0])
            mosblock = []
            mosblock.append('   '+str(nalpha+idx)+' '+words[1]+' '+words[2]+' '+words[3]+'\n')
            line = fortfile.readline()
            while 'MO' not in line:
                mosblock.append(line)
                line = fortfile.readline()
            nb += 1
            betamos.append(mosblock)
    fortfile.close()
    print("fort7 file read")

    # output redirection
    sys_stdout_orig = sys.stdout
    sys.stdout = open(wffilename,'w')
    print('$general')
    print("title='"+title+"'")
    print('charge='+str(charge)+', spin='+str(mult)+', geom=angstrom')
    print('evfmt=gau, basis='+basis+', jastrow=none')
    if charge!=0:
      print('atomic_charges')
    print('$end')
    print('$geom')
    print(len(geometry))
    pcharge = [0]*len(geometry)    
    if charge!=0:
      eneg = [[0.0,0]]*len(geometry)    
      for i in range(len(geometry)):
        eneg[i] = [electronegativity[pse.index(geometry[i][0])],i]
      if charge >0 :
        #positive overall charge, remove electrons
        tcharge = charge
        while tcharge != 0:
          tidx = sorted(eneg, key=get1st)[0][1]
          if (eneg[tidx][0] != eneg[tidx][0]):
             # lowest EN is 'nan' -> all EN are nan
             raise NameError('Wrong charge given.')
          pcharge[tidx] = pcharge[tidx] + 1
          eneg[tidx][0] = eneg[tidx][0]+0.5
          if (pse.index(geometry[i][0])-pcharge[tidx] <= 0):
            eneg[tidx][0] = float('nan') 
          tcharge = tcharge - 1
      else:
        #negative overall charge, add electrons        
        tcharge = charge
        while tcharge != 0:
          tidx =  sorted(eneg, key=get1st, reverse=True)[0][1]
          pcharge[tidx] = pcharge[tidx] - 1
          eneg[tidx][0] = eneg[tidx][0]-0.5
          tcharge = tcharge + 1
    idx = -1                     
    for line in geometry:
      idx = idx+1
      if charge!=0:
        if basis=='diff':
          if line[0] in ecpatoms:
            print(line[0].rjust(2),'   ',line[2],line[3],line[4],str(pcharge[idx]),' ',ecp)
          else:
            print(line[0].rjust(2),'   ',line[2],line[3],line[4],str(pcharge[idx]),' ',aebasis)
        else:
          print(line[0].rjust(2),'   ',line[2],line[3],line[4],str(pcharge[idx]))
      elif basis=='diff':
        if line[0] in ecpatoms:
          print(line[0].rjust(2),'   ',line[2],line[3],line[4],' ',ecp)
        else:
          print(line[0].rjust(2),'   ',line[2],line[3],line[4],' ',aebasis)
      else:
        print(line[0].rjust(2),'   ',line[2],line[3],line[4])
    print('$end')
    print('$mos')
    if method[0] == 'U':
        print(nalphaorb+nbetaorb)
    else:
        print(nalphaorb)
    print(' ')
    for e in alphamos:
        for line in e:
            print(line[:-1])
    if method[0] == 'U':
        for e in betamos:
            for line in e:
                print(line[:-1])
    print('$end')
    print('$dets')
    #print('1')
    #s = '  1.0  '
    #for a in range(1,nalpha+1):
    #    s += ' '+str(a)
    #s += ' '
    #if method[0] == 'U':
    #    for b in range(1,nbeta+1):
    #        s += ' '+str(nalpha+b)
    #else:
    #    for a in range(1,nbeta+1):
    #        s += ' '+str(a)
    #print(s)
    if findDets:
        used_dets = sorted(
            [(ev, d) for ev, d in zip(eigenvectors, det) if abs(ev) >= coeffcutoff],
            key=lambda a:-abs(a[0])
        )
        print(len(used_dets))
        for ev, det in used_dets:
            print("%s %s" % (ev, det))
    else:
        if method[0] == 'U':
            print('single unrestricted')
        else:
            print('single restricted')
    print('$end')
    sys.stdout.close()
    sys.stdout = sys_stdout_orig
    print(' amolqc wf file '+wffilename+' has been created')


def readGeometry(logfile, orient, pse):
    logfile.seek(0)
    norients = 0
    line = ' '
    while line:
        if orient in line:
            norients += 1
        line = logfile.readline()
    if norients == 0: raise EOFError("geometry not found in log file")
    logfile.seek(0)
    while True:
        line = logfile.readline()
        if orient in line:
            norients -= 1
            if norients == 0:
                for i in range(5):
                    line = logfile.readline()
                geometry = []
                while '----' not in line:
                    words = line.split()
                    geometry.append([pse[int(words[1])], int(words[1]), words[3], words[4], words[5]])
                    line = logfile.readline()
                break
    return geometry


#
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('''    gaussian2wf.py creates an amolqc wf file from
    gaussian output file and the fort.7 file
    usage:
    gaussian2wf.py outfile -b basisname [-m method ][-t title] [-ns]
    outfile:  without extension
    method: 'R' or 'U' (or full method name), default is 'R'
    title: title to include in wf file, default is 'no title'
    -ns: NoSymm calculation, read Input orientation
    -c: Cutoff for determinant coefficients
    -nd: Don't extract determinants
    ''')
        sys.exit(1)

    fname = sys.argv[1]
    args = sys.argv[2:]
    basis = 'XXX'
    title = 'no title'
    method = 'R'
    coeffcutoff = 0.0000001
    stdorient = True
    dets = True
    while len(args)>0:
        if args[0]=='-t':
            title = args[1]
            args = args[2:]
        elif args[0]=='-b':
            basis = args[1]
            args = args[2:]
        elif args[0]=='-m':
            method = args[1]
            args = args[2:]
        elif args[0]=='-ns':
            stdorient = False
            args = args[1:]
        elif args[0]=='-nd':
            dets = False
            args = args[1:]
        elif args[0] == '-c':
            coefcutoff = float(args[1])
            args = args[2:]
        else:
            print('unrecognized option %s' % args[0])
            sys.exit(1)

    gaussian2wf(fname,basis,method=method,title=title,stdorient=stdorient,coeffcutoff=coeffcutoff,
        dets=dets)
