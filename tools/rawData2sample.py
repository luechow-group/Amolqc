#!/usr/bin/env python
#
""""
 Copyright (C) 2020 Leonard Reuter

 SPDX-License-Identifier: GPL-3.0-or-later
"""

import sys
if sys.version_info[0] < 3:
    sys.exit('This script requires Python 3')

import numpy as np
import os

i4 = np.int32
r8 = np.float64


def fortranify(tpllist, everywhere=True):
    # by default before and after every entry, fortran saves the length of the entry
    # as a 4 byte integer. 'fortranify' converts a tpllist to a dtype accordingly.

    ftpllist = []
    if everywhere:
        for i, element in enumerate(tpllist, 1):
            ftpllist.append((f'{i}A', i4))
            ftpllist.append(element)
            ftpllist.append((f'{i}E', i4))
    else:
        ftpllist.append((f'A', i4))
        for i, element in enumerate(tpllist, 1):
            ftpllist.append(element)
        ftpllist.append((f'E', i4))
    return np.dtype(ftpllist)


def read_rawdata(base_name, mode, verbose=True, atoms=False, nalpha=False):
    vector = ''
    if mode == 'maxima':
        vector = 'maximum positions vector'
    elif mode == 'sample':
        vector = 'sample positions vector'
    else:
        sys.exit("Error: sample has to be 'maxima' or 'sample'.")

    i = 0
    file_name = f"{base_name}-{i:02d}.bin"
    if not os.path.isfile(file_name):
        sys.exit(file_name + " does not exist")

    # determining number of atoms
    N_dtype = fortranify([('number of atoms', i4)])
    N = np.fromfile(file_name, dtype=N_dtype, count=1)['number of atoms'][0]

    # determining number of electrons
    n_dtype = fortranify([
        ('number of atoms', i4),
        ('atomic numbers vector', i4, (N)),
        ('atom positions vector', r8, (3*N)),
        ('number of electrons', i4)])
    n = np.fromfile(file_name, dtype=n_dtype, count=1)['number of electrons'][0]

    # defining dtypes
    header_tpllist = [
        ('number of atoms', i4),
        ('atomic numbers vector', i4, (N)),
        ('atom positions vector', r8, (3*N)),
        ('number of electrons', i4),
        ('number of alpha_electrons', i4)]

    header_dtype = fortranify(header_tpllist)

    body_tpllist = [
        ('sample positions vector', r8, (3*n)),
        ('kinetic energies vector', r8, (n)),
        ('maximum positions vector', r8, (3*n)),
        ('-ln|Psi|^2 value', r8)]

    body_dtype = fortranify(body_tpllist, everywhere=False)

    # reading sample as maxima or sample from rawdata file
    sample = np.fromfile(file_name, dtype=body_dtype, count=-1, offset=header_dtype.itemsize)[vector]
    if verbose:
        print(f'read {file_name}')
    i = 1
    file_name = f"{base_name}-{i:02d}.bin"
    while os.path.isfile(file_name):
        sample = np.append(sample, np.fromfile(file_name, dtype=body_dtype, count=-1)[vector], axis=0)
        if verbose:
            print(f'read {file_name}')
        i += 1
        file_name = f"{base_name}-{i:02d}.bin"

    if verbose:
        print(f'total sample size: {len(sample)}')

    results = [sample, n]
    if atoms:
        header = np.fromfile(f"{base_name}-00.bin", dtype=header_dtype, count=1)
        results += [header['atomic numbers vector'][0], header['atom positions vector'][0]]
    if nalpha:
        header = np.fromfile(f"{base_name}-00.bin", dtype=header_dtype, count=1)
        results += [header['number of alpha_electrons'][0]]

    return tuple(results)


def write_sample(sample, n, sample_filename='sample.pos', files_number=1):
    # writing sample to binary (again with fortranic length markers)
    header_length = np.array([8], dtype=i4)
    sample_length = np.array([3*8*n], dtype=i4)
    n_array = np.array([n], dtype=i4)
    tsize_array = np.array([len(sample)], dtype=i4)
    with open(sample_filename, 'wb') as outfile:
        outfile.write(header_length.tobytes())
        for array in [n_array, tsize_array]:
            outfile.write(array.tobytes())
        outfile.write(header_length.tobytes())
        for i, walker in enumerate(sample):
            outfile.write(sample_length.tobytes())

            # each walker has to be rearranged
            # from xyzxyzxyz
            # to xxxyyyzzz
            for j in range(3):
                mask = np.array(range(3*n)) % 3 == j
                outfile.write(walker[mask].tobytes())
            outfile.write(sample_length.tobytes())
            # if (i+1) % (len(sample) // 100) == 0:
            #    print(f'wrote {(i+1)//(len(sample) // 100)}% to {sample_filename}')
    print(f'wrote {sample_filename}')


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(
            '''
            
                Creates sample binary from rawData binary and prints the sample size.
                
                Usage: rawData2sample.py rawBaseName vector [file_number]
            
                rawBaseName:       base name of the rawData binaries (eg. 'raw' for 'raw-00.bin')
                vector:            either 'maxima' or 'sample', depending on what should be read from the rawData file
                file_number:       sample is split among this number of files for writing, default=1 (ie. no splitting)
            
            ''')

    base_name_ = sys.argv[1]
    mode_ = sys.argv[2]

    sample_basename = f'{mode_}'

    files_number = 1
    if len(sys.argv) > 3:
        files_number = int(sys.argv[3])

    sample_, n_ = read_rawdata(base_name_, mode_)

    nf = len(sample_) // files_number
    print(f'sample size per file: {nf}')
    for k in range(files_number):
        sample_filename_ = f'{sample_basename}-{k+1:02d}.pos'
        write_sample(sample_[k*nf:(k+1)*nf], n_, sample_filename_)
