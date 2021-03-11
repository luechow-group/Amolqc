"""
 Copyright (C) 2017-2020 Leonard Reuter
 Copyright (C) 2020 Arne Luechow

SPDX-License-Identifier: GPL-3.0-or-later
"""

from sys import exit, stdout
import itertools as it

from .csf import Csf
from .orbitalRotation import OrbitalRotation
from .data import multiplication_table, periodic_table


class WaveFunction:
    def __init__(self, input_name, orbital_format, basis, charge, multiplicity, atoms, orbitals,
                 csfs, jastrow, symmetry=None, symmetry_list=None):
        self.title = input_name  # string
        self.orbital_format = orbital_format  # string
        self.basis = basis  # string
        self.charge = charge  # integer
        self.multiplicity = multiplicity  # integer
        self.atoms = atoms  # list of atom objects
        self.orbitals = orbitals  # list of orbital objects
        self.csfs = csfs  # list of csf objects
        self.symmetry = symmetry  # string
        self.jastrow = jastrow
        self.symmetry_list = symmetry_list
        self.separated_electrons = 0

    def write(self):
        self.title = self.title.split('/')[-1]
        wf_filename = self.title + '.wf'
        out = open(wf_filename, 'w')
        self.write_general_section(out)

        # geometry
        self.write_geometry_section(out)

        # jastrow
        if self.jastrow.type != 'none':
            self.write_jastrow_section(out)

        # basis
        if self.basis == 'gaussian':
            self.write_basis_section(out)

        # orbitals
        self.write_orbital_section(out)

        if len(self.csfs) != 1:
            # csfs
            self.write_csfs(out)
        else:
            self.write_dets(out)

        out.close()

    def write_jastrow_section(self, out):
        out.write('$jastrow' + '\n')
        self.jastrow.write(out)
        out.write('$end' + '\n')

    def write_general_section(self, out):
        out.write('$general' + '\n'
                  + ' title=' + self.title + '\n'
                  + ' evfmt=' + self.orbital_format + ', basis=' + self.basis + ', jastrow=' + self.jastrow.type + '\n'
                  + ' charge= ' + str(self.charge) + ', spin= ' + str(self.multiplicity))
        if self.atomic_charges():
            out.write(', atomic_charges')
        if self.separated_electrons != 0:
            out.write(f',\n separated_electrons={self.separated_electrons}')
        out.write('\n'
                  + '$end' + '\n')

    def write_geometry_section(self, out):
        out.write('$geom' + '\n'
                  + '  ' + str(len(self.atoms)) + '\n')
        for atom in self.atoms:
            atom.write(out, self.atomic_charges())
        out.write('$end' + '\n')

    def write_basis_section(self, out):
        out.write('$basis' + '\n')
        for a in self.atoms:
            out.write(periodic_table[a.atomic_number] + '\n')
            for bf in a.basis_functions:
                bf.write_amolqc(out)
            out.write(' ****' + '\n')
        out.write('$end' + '\n')

    def write_orbital_section(self, out):
        out.write('$mos' + '\n'
                  + '  ' + str(len(self.orbitals)) + '\n\n')
        i = 0
        for orb in self.orbitals:
            orb.write(out, i, self.orbital_format)
            i += 1
        out.write('$end' + '\n')

    def write_csfs(self, out):
        out.write('$csfs' + '\n'
                  + str(len(self.csfs)).rjust(7) + '\n')
        for csf in self.csfs:
            csf.write(out)
        out.write('$end' + '\n')

    def write_dets(self, out):
        if len(self.csfs[0].determinants) == 1:
            self.csfs[0].determinants[0].coefficient = 1.0
        out.write('$dets' + '\n'
                  + str(len(self.csfs[0].determinants)).rjust(7) + '\n')
        for det in self.csfs[0].determinants:
            det.write(out)
        out.write('$end' + '\n')

    def print_orbital_symmetries(self):
        i = 0
        for orb in self.orbitals:
            print(str(i+1)+' '+str(orb.symmetry))
            i += 1

    def punch_vec(self):
        vec = open(self.title + '.vec', 'w')
        vec.write(' $VEC' + '\n')
        i = 0
        for orb in self.orbitals:
            orb.write(vec, i, 'gms')
            i += 1
        vec.write(' $END' + '\n')
        vec.close()

    def sort_dets(self):
        if len(self.csfs) == 1 and len(self.csfs[0].determinants) == 1:
            print('Warning: sort_dets does nothing for single determinant wave functions.')

        for csf in self.csfs:
            for det in csf.determinants:
                det.sort()

    def convert_to_det(self):
        if len(self.csfs) == 1:
            print('Warning: convert does nothing if the number of csfs is already 1.')

        self.sort_dets()
        old_csfs = self.csfs
        self.csfs = []
        new_csf = Csf()
        for csf in old_csfs:
            for det in csf.determinants:
                new_det = det
                new_det.coefficient *= csf.coefficient
                if new_det in new_csf.determinants:
                    new_csf.determinants[new_csf.determinants.index(new_det)].coefficient += new_det.coefficient
                else:
                    new_csf.determinants.append(new_det)
        self.csfs.append(new_csf)

    def sort_ci(self):
        if len(self.csfs) == 1 and len(self.csfs[0].determinants) == 1:
            print('Warning: sort_ci does nothing for single determinant wave functions.')

        for csf in self.csfs:
            csf.determinants.sort(key=lambda determinant: -1.0 * abs(determinant.coefficient))
        self.csfs.sort(key=lambda csf: -1.0 * abs(csf.coefficient))

    def alt_sort_ci(self):
        if len(self.csfs) == 1 and len(self.csfs[0].determinants) == 1:
            print('Warning: alt_sort_ci does nothing for single determinant wave functions.')

        for csf in self.csfs:
            csf.determinants.sort(key=lambda determinant: determinant.occupation_integer())
        # self.csfs.sort(key=lambda csf: -1.0 * abs(csf.coefficient))
        self.csfs.sort(key=lambda csf: csf.determinants[0].occupation_integer())

    def cut_ci(self, threshold):
        self.sort_ci()
        i = 0
        if len(self.csfs) > 1:
            while abs(self.csfs[i].coefficient) > threshold:
                i += 1
                if i == len(self.csfs):
                    print('Warning: no CI coeffs are below the given threshold, the wf is unchanged.')
                    break
            for _ in range(len(self.csfs) - i):
                del self.csfs[i]
        else:
            print('Warning: dets are cut, since only one csf is found. The wf is no longer spin symmetric.')
            while abs(self.csfs[0].determinants[i].coefficient) > threshold:
                i += 1
                if i == len(self.csfs[0].determinants):
                    print('Warning: no CI coeffs are below the given threshold, the wf is unchanged.')
                    break
            for _ in range(len(self.csfs[0].determinants) - i):
                del self.csfs[0].determinants[i]

    def symm_combine(self):
        if len(self.csfs) == 1:
            print('Warning: symm_combine does nothing if the number of csfs is already 1.')

        question_string = "How many distinct symmetry combinations should be performed? "
        input_string = input(question_string + "(Enter h for help.) ")
        if input_string == 'h':
            input_string = input("\nIn one distinct symmetry combination, several orbital exchanges can be checked simultaneously.\n"
                                 "This is for instance useful for pi-symmetry, since pi_x -> pi_y and pi*_x -> pi*_y should only\n"
                                 "be checked at the same time!\n"
                                 "Several distinct symmetry combinations may occur, if pi and delta symmetry should both be\n"
                                 "considered since they are independent from one another.\n\n"
                                 + question_string)
        count = int(input_string)
        for ii in range(count):
            self.sort_dets()
            counter = 0
            orbital_groups = self.automatic_orbital_groups()

            list_for_permutation = range(len(orbital_groups[0]))
            permutations = list(it.permutations(list_for_permutation))
            del permutations[0]
            i = 0
            while i < len(self.csfs):
                j = i + 1
                for k in range(len(permutations)):
                    permuted_occupation = list(self.csfs[i].occupation)
                    if self.symmetry != 'VB':
                        for l in range(len(orbital_groups[0])):
                            if permutations[k][l] != l:
                                for m in range(len(orbital_groups)):
                                    if len(orbital_groups[m]) <= permutations[k][l]:
                                        exit('Error: symm based csf combination failed, check input groups or try coeff_combine.')
                                    permuted_occupation[orbital_groups[m][l]] = self.csfs[i].occupation[orbital_groups[m][permutations[k][l]]]
                    else:
                        # TODO
                        raise NotImplementedError('VB symmetry combination is not implemented')

                    # Build alternative permuted occupation exchanging "/" and "\"
                    alt_permuted_occupation = list(permuted_occupation)
                    for l in range(len(alt_permuted_occupation)):
                        if alt_permuted_occupation[l] == '/':
                            alt_permuted_occupation[l] = '\\'
                        elif alt_permuted_occupation[l] == '\\':
                            alt_permuted_occupation[l] = '/'

                    while j < len(self.csfs):
                        if (permuted_occupation == self.csfs[j].occupation
                            or alt_permuted_occupation == self.csfs[j].occupation) \
                                and len(self.csfs[i].determinants) % len(self.csfs[j].determinants) == 0 \
                                and len(self.csfs[i].determinants) // len(self.csfs[j].determinants) < len(orbital_groups[0]):
                            if not abs(abs(self.csfs[i].coefficient) - abs(self.csfs[j].coefficient)) < 1E-6:
                                print('Warning: based on the symmetry, two csfs with > 1E-6 differences in coefficient'
                                      ' are to be merged.')
                                combine = self.ask_csf_combination(i, j)
                                if combine:
                                    self.merge_two_csfs(i, j)
                                    counter += 1
                                else:
                                    j += 1
                            else:
                                self.merge_two_csfs(i, j)
                                counter += 1
                        else:
                            j += 1
                i += 1

            print(str(counter)+" csfs merged based on symmetry.")
            self.check_csfs_coeffs()

            for csf in self.csfs:
                csf.normalize()

    def automatic_orbital_groups(self):
        print('''
        Give the the indices of degenerate irreps,
        separated with commas. When finished, press enter once again.''')
        print(self.symmetry)
        for i in range(len(self.symmetry_list)):
            print(str(i + 1) + ': ' + self.symmetry_list[i])
        group_list = []
        while True:
            input_list = input('Group of degenerate irreps: ').split(',')
            if len(input_list) != 1:
                group_list.append(input_list)
            else:
                break
        orbital_groups = []
        for i in range(len(group_list)):
            for j in range(len(self.orbitals)):
                if self.orbitals[j].symmetry.split('.')[1] == group_list[i][0]:
                    orbital_group = [j]
                    for k in range(len(group_list[i]) - 1):
                        for l in range(len(self.orbitals)):
                            if self.orbitals[l].symmetry.split('.')[1] == group_list[i][k + 1] \
                            and self.orbitals[l].symmetry.split('.')[0] == self.orbitals[j].symmetry.split('.')[0]:
                                orbital_group.append(l)
                    orbital_groups.append(orbital_group)
        return orbital_groups

    def get_virtual_orbitals(self):
        virtual_orbitals = []
        for i in range(len(self.orbitals)):
            index = i + 1
            empty = True
            for csf in self.csfs:
                for determinant in csf.determinants:
                    # print(determinant.orbital_list)
                    if not (determinant.orbital_list.count(index) == 0
                            and determinant.orbital_list.count(-1 * index) == 0):
                        empty = False
                        break
                if not empty:
                    break
            if empty:
                virtual_orbitals.append(index)

        return virtual_orbitals

    def get_inactive_orbitals(self):
        inactive_orbitals = []
        for i in range(len(self.orbitals)):
            index = i + 1
            occupied = True
            for csf in self.csfs:
                for determinant in csf.determinants:
                    # print(determinant.orbital_list)
                    if not (determinant.orbital_list.count(index) == 1
                            and determinant.orbital_list.count(-1 * index) == 1):
                        occupied = False
                        break
                if not occupied:
                    break
            if occupied:
                inactive_orbitals.append(index)

        return inactive_orbitals

    def get_active_orbitals(self):
        active_orbitals = []
        inactive_orbitals = self.get_inactive_orbitals()
        virtual_orbitals = self.get_virtual_orbitals()

        for i in range(len(self.orbitals)):
            index = i + 1
            if index not in inactive_orbitals + virtual_orbitals:
                active_orbitals.append(index)
        return active_orbitals

    def punch_moopt(self):
        file = open(self.title + '.moopt', 'w')

        # getting relevant orbital indices
        inactive_orbitals = self.get_inactive_orbitals()
        virtual_orbitals = self.get_virtual_orbitals()
        active_orbitals = self.get_active_orbitals()

        # get orbital_rotations
        orbital_rotations = []
        orbital_rotations += self.get_orbital_rotations(inactive_orbitals, active_orbitals)
        orbital_rotations += self.get_orbital_rotations(inactive_orbitals, virtual_orbitals)
        orbital_rotations += self.get_orbital_rotations(active_orbitals, virtual_orbitals)

        query = input('Should non-required virtual orbitals be deleted? [y/N]')
        if query.lower() in ['y', 'yes']:
            # identifying orbitals, that are in no orbital_rotation
            orbital_mask = [True]*len(self.orbitals)
            for i in range(len(self.orbitals)):
                for orbital_rotation in orbital_rotations:
                    if (i+1 in orbital_rotation.orbital_group[0]) or (i+1 in orbital_rotation.orbital_group[1]):
                        break
                else:  # no break
                    orbital_mask[i] = False
            index_mask = [index for index, required in zip(range(len(self.orbitals)), orbital_mask) if required]

            # removing non-required orbitals
            self.orbitals = [orbital for orbital, required in zip(self.orbitals, orbital_mask) if required]

            # adapting orbital_rotations
            for orbital_rotation in orbital_rotations:
                for i in [0, 1]:
                    orbital_rotation.orbital_group[i] = [index_mask.index(orb-1)+1
                                                         for orb in orbital_rotation.orbital_group[i]]
        elif query.lower() in ['n', 'no', '']:
            pass
        else:
            exit('pleaser answer y or n')

        file.write('orbital_rotation_list=\n')
        file.write(str(len(orbital_rotations)) + '\n')

        for excitation in orbital_rotations:
            excitation.write(file)

        orbital_groups = self.automatic_orbital_groups()

        file.write('\n')
        file.write('mo_symmetrise_list=\n')
        file.write(str(len(orbital_groups)) + '\n')

        for orbital_group in orbital_groups:
            file.write(str(len(orbital_group)))
            for index in orbital_group:
                file.write(' ' + str(index + 1))
            file.write('\n')

        file.close()

    def coeff_combine(self):
        print('Warning: coeff_combine should be used very carefully. Please check every combination.')
        if len(self.csfs) == 1:
            print('Warning: coeff_combine does nothing if the number of csfs is already 1.')

        i = 0
        change = False
        while i < len(self.csfs):
            j = i + 1
            combined = False
            while j < len(self.csfs):
                if abs(abs(self.csfs[i].coefficient) - abs(self.csfs[j].coefficient)) < 1E-6 \
                        and len(self.csfs[i].determinants) % len(self.csfs[j].determinants) == 0 \
                        and (len(self.atoms) == 1 or len(self.csfs[i].determinants)
                             // len(self.csfs[j].determinants) == 1):
                    if combined:
                        print("WARNING: CSF A is already a combined CSF!!!")
                    combine = self.ask_csf_combination(i, j)
                    if combine:
                        self.merge_two_csfs(i, j)
                        combined = True
                        change = True
                    else:
                        j += 1
                else:
                    j += 1
            i += 1
        if not change:
            print("Warning: No CSFs have been combined.")

        for csf in self.csfs:
            csf.normalize()

    def force_combine(self):
        self.sort_dets()
        # print("Warning: CSFs are combined based on coefficients without checking.")
        i = 0
        # change = False
        while i < len(self.csfs):
            j = i + 1
            while j < len(self.csfs):
                if abs(abs(self.csfs[i].coefficient) - abs(self.csfs[j].coefficient)) < 1E-6:
                    self.merge_two_csfs(i, j)
                    # change = True
                else:
                    j += 1
            i += 1
        # if not change:
        #    print("Warning: No CSFs have been combined.")
        for csf in self.csfs:
            csf.normalize()

    def check_csfs_coeffs(self):
        counter = 0
        i = 0
        while i < len(self.csfs):
            j = i + 1
            while j < len(self.csfs):
                if abs(abs(self.csfs[i].coefficient) - abs(self.csfs[j].coefficient)) < 1E-6:
                    counter += 1
                j += 1
            i += 1
        print("Warning: based on the coeff difference (< 1E-6), "+str(counter)+" more combinations are thinkable.\n"
              + "You may want to carefully check them with 'wfgen.py amolqc *.wf coeff_combine write'.")

    def ask_csf_combination(self, i, j):
        combine = False
        print("CSF A")
        self.csfs[i].write(stdout)
        print("CSF B")
        self.csfs[j].write(stdout)

        test = input("Should these csfs be combined? (Y/n): ").lower()
        if test == 'yes' or test == 'y' or test == '':
            combine = True
            print('<<< CSFS COMBINED >>>')
        elif test != 'no' and test != 'n':
            exit("Error: user input '" + test + "' is not valid.")
        print('---------------------')

        return combine

    def merge_two_csfs(self, i, j):
        if abs(self.csfs[i].coefficient - self.csfs[j].coefficient) \
                < abs(self.csfs[i].coefficient + self.csfs[j].coefficient):
            sign = 1
        else:
            sign = -1
        for k in range(len(self.csfs[j].determinants)):
            if self.csfs[j].determinants[k] in self.csfs[i].determinants:
                index = self.csfs[i].determinants.index(self.csfs[j].determinants[k])
                self.csfs[i].determinants[index].coefficient += sign * self.csfs[j].determinants[k].coefficient
                if abs(self.csfs[i].determinants[index].coefficient) < 1E-7:
                    del self.csfs[i].determinants[index]
            else:
                determinant = self.csfs[j].determinants[k]
                determinant.coefficient *= sign
                self.csfs[i].determinants.append(determinant)
        del self.csfs[j]

    def rename(self, name):
        self.title = name+''  # to check variable type

    def set_charges(self):
        if len(self.atoms) > 1:
            print('The charge is ' + str(self.charge) + '.')
            for i in range(len(self.atoms)):
                stdout.write(str(i+1).rjust(2) + ': ')
                self.atoms[i].write(stdout, self.atomic_charges())
            index_string = input('Give the line of the atom which bears the charge: ')
            index = 0
            try:
                index = int(index_string)
            except ValueError:
                exit('Error: the line number must be an integer!')
            for i in range(len(self.atoms)):
                if i == index - 1:
                    self.atoms[i].charge = self.charge
                else:
                    self.atoms[i].charge = 0
        else:
            self.atoms[0].charge = self.charge

    def atomic_charges(self):
        for atom in self.atoms:
            if atom.charge != 0:
                return True
        else:
            return False

    def check_normalization(self):
        normalization = self.calculate_normalization()
        print('Integral = '+str(normalization) + '. (Should be 1.0)')

    def check_coeff_sum(self):
        coeff_sum = self.calculate_coeff_sum()
        print('Sum of abs coeffs = '+str(coeff_sum) + '. (Should be > 1.0)')

    def check_excitation(self):
        excitation = self.calculate_excitation()
        coeff_sum = self.calculate_coeff_sum()
        print('Excitation = '+str(excitation) + '. (of '+str(coeff_sum)+')')

    def calculate_normalization(self):
        coeff_sum = 0.0
        for csf in self.csfs:
            csf_sum = 0.0
            for det in csf.determinants:
                csf_sum += det.coefficient ** 2
            coeff_sum += csf_sum * (csf.coefficient ** 2)
        return coeff_sum

    def calculate_coeff_sum(self):
        coeff_sum = 0.0
        for csf in self.csfs:
            csf_sum = 0.0
            for det in csf.determinants:
                csf_sum += abs(det.coefficient)
            coeff_sum += csf_sum * (abs(csf.coefficient))
        return coeff_sum

    def calculate_excitation(self):
        orbitals = input('Give orbital indices separated by commas: ').split(',')
        excitation_sum = 0.0
        for csf in self.csfs:
            csf_sum = 0.0
            for det in csf.determinants:
                excitation = False
                for orb in orbitals:
                    if (det.orbital_list.count(int(orb)) + det.orbital_list.count(-int(orb))) != 2:
                        excitation = True
                if excitation:
                    csf_sum += det.coefficient ** 2
            excitation_sum += csf_sum * (csf.coefficient ** 2)
        return excitation_sum

    def normalize(self):
        for csf in self.csfs:
            csf.normalize()
        coeff_sum = 0.0
        for csf in self.csfs:
            coeff_sum += csf.coefficient ** 2
        factor = coeff_sum ** 0.5
        for csf in self.csfs:
            csf.coefficient /= factor

    def cut_orbs(self):
        max_orb = 0
        for csf in self.csfs:
            for det in csf.determinants:
                for orb in det.orbital_list:
                    max_orb = max(abs(orb), max_orb)
        initial_orbs = len(self.orbitals)
        for i in range(initial_orbs - max_orb):
            del self.orbitals[max_orb]

    def max_dets(self):
        max_det = 0
        for csf in self.csfs:
            max_det = max(max_det, len(csf.determinants))
        print("Maximal number of determinants per csf is " + str(max_det) + ".")

    def add_symmetry(self):
        count = int(input("How many new symmetries should be added? "))
        self.symmetry += "+"
        for i in range(count):
            symm_name = input('Give the name for the symmetry '+str(i + 1)+': ')
            self.symmetry += symm_name
            if i < count - 1:
                self.symmetry += "/"
            self.symmetry_list.append(symm_name)

            orbital_list = input('Give orbital indices, that should be added to '
                                 + symm_name + ', seperated by comma: ').split(",")

            counter = 0
            for j in range(len(self.orbitals)):
                if str(j + 1) in orbital_list:
                    counter += 1
                    self.orbitals[j].symmetry = str(counter) + '.' + str(len(self.symmetry_list))
            # print(self.symmetry_list)
            print(str(counter) + ' orbitals have been added to ' + self.symmetry_list[-1])

    def get_orbital_rotations(self, group1_orbitals, group2_orbitals):
        orbital_rotations = []
        for i, irrep in enumerate(self.symmetry_list):
            group1 = []
            group2 = []
            for index in group1_orbitals:
                if self.orbitals[index - 1].symmetry.split('.')[1] == str(i+1):
                    group1.append(index)
            for index in group2_orbitals:
                if self.orbitals[index - 1].symmetry.split('.')[1] == str(i+1):
                    group2.append(index)
            if len(group1) != 0 and len(group2) != 0:
                orbital_rotation = OrbitalRotation()
                orbital_rotation.orbital_group.append(group1)
                orbital_rotation.orbital_group.append(group2)
                orbital_rotations.append(orbital_rotation)
        return orbital_rotations

    def print_csf(self, index):
        assert(index > 0)
        assert(index <= len(self.csfs))
        self.csfs[index - 1].write(stdout)

    def get_csfs_sharing_determinants(self, verbose=True):
        if verbose:
            print('Pairs of csfs sharing determinants:')
        sharing = False
        for i in range(len(self.csfs)):
            for j in range(i+1, len(self.csfs)):
                for determinant in self.csfs[i].determinants:
                    if determinant in self.csfs[j].determinants:
                        if verbose:
                            print(i, j)
                        sharing = True
                        break

        if verbose and not sharing:
            print('None')

        return sharing

    def replace_orbital_indices(self):
        replacements = 0
        affectedDets = []
        oldIndices = input('Give orbitals to replace in dets, seperated by spaces: ').split()
        newIndices = input('Give new indices in the respective oder: ').split()
        for indices in zip(oldIndices, newIndices):
            oldIndex = int(indices[0])
            newIndex = int(indices[1])
            for i, csf in enumerate(self.csfs):
                for j, det in enumerate(csf.determinants):
                    detReplaced = False
                    for k, index in enumerate(det.orbital_list):
                        if abs(index) == oldIndex:
                            if index == oldIndex:
                                self.csfs[i].determinants[j].orbital_list[k] = newIndex
                            else:
                                self.csfs[i].determinants[j].orbital_list[k] = -newIndex
                            replacements += 1
                            detReplaced = True
                    if detReplaced:
                        detTuple = (i,j)
                        if detTuple not in affectedDets:
                            affectedDets.append(detTuple)
        print(f'Replaced {replacements} orbital indices in {len(affectedDets)} determinants.')

    def reorder_atoms(self):
        atom_types = set([atom.atomic_number for atom in self.atoms])

        bf_numbers = {}
        print("Reordering atoms:")
        print(f"The total number of basis functions is {len(self.orbitals[0].coefficients)}.")
        for atom_type in atom_types:
            bf_numbers[atom_type] = int(input(f"Number of basis functions for element {periodic_table[atom_type]}: "))

        assert(sum([bf_numbers[atom_type] for atom_type in [atom.atomic_number for atom in self.atoms]])
               == len(self.orbitals[0].coefficients)), "Error in number of basis functions."

        input_list = input('Give the new atom order, comma seperated: ').split(',')
        reorder_list = [int(i)-1 for i in input_list]

        # make sure, no atom is twice in reorder list:
        assert(len(set(reorder_list)) == len(reorder_list)), "At least one atom appeared twice in the new atom order."
        assert(set(reorder_list).issubset(set(range(len(self.atoms))))), "At least one index is invalid."

        # if reorder_list is shorter than len(self.atoms), all other atoms are implicitly not changed
        if len(reorder_list) != len(self.atoms):
            print("Not all atom indices appear in the new atom order. The non-appearing atoms are implicitly assumed"
                  " to be keeping their index.")
            for index in range(len(self.atoms)):
                if index not in reorder_list:
                    reorder_list.insert(index, index)

        # reordering orbitals
        for i, orbital in enumerate(self.orbitals):
            coefficients = []
            for r in reorder_list:
                start_index = 0
                for index in range(r):
                    start_index += bf_numbers[self.atoms[index].atomic_number]
                end_index = start_index + bf_numbers[self.atoms[r].atomic_number]
                coefficients += orbital.coefficients[start_index:end_index]
            self.orbitals[i].coefficients = coefficients

        # reordering atoms
        self.atoms = [self.atoms[r] for r in reorder_list]

    def separate_electrons(self):
        print("Separating electrons:")
        self.separated_electrons = int(input("Give the number of separated electrons: "))
        assert(self.separated_electrons % 2 == 0), "With the current implementation, separated electrons has to be even"

        for i, csf in enumerate(self.csfs):
            for j, det in enumerate(csf.determinants):
                self.csfs[i].determinants[j].orbital_list[:] = (orb for orb in det.orbital_list
                                                                if abs(orb) > self.separated_electrons // 2)

        atom_types = set([atom.atomic_number for atom in self.atoms])
        charges = {}
        for atom_type in atom_types:
            charges[atom_type] = int(input(f"Charge (for sampling) for element {periodic_table[atom_type]}: "))

        for i, atom in enumerate(self.atoms):
            self.atoms[i].charge = charges[atom.atomic_number]
