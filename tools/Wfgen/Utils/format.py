""""
 Copyright (C) 2017-2019 Leonard Reuter

SPDX-License-Identifier: GPL-3.0-or-later
"""

from sys import exit


def format_coordinate(coord):
    formatted_coordinate = '{: 06.5f}'.format(coord)
    return formatted_coordinate


def format_det_coeff(coeff):
    formatted_coeff = '{: 08.7f}'.format(coeff)
    return formatted_coeff


def format_orbital_coefficient(coefficient, orbital_format):
    if orbital_format == 'gms':
        formatted_coefficient = '{: 09.8E}'.format(coefficient)
    else:
        if abs(float(coefficient)) > 0:
            formatted_coefficient = '{: 08.7E}'.format(coefficient)
            temp_list = list(formatted_coefficient)
            temp_list.insert(3, temp_list[1])
            temp_list[1] = '0'
            temp_int = int(''.join([temp_list[12], temp_list[13], temp_list[14]]))
            for i in range(3):
                del temp_list[12]
            formatted_coefficient = ''.join(temp_list)
            temp_int += 1
            formatted_coefficient += '{0:+03d}'.format(temp_int)
        else:
            formatted_coefficient = '{: 09.8E}'.format(coefficient)
        formatted_coefficient = formatted_coefficient.replace("E", "D")
    return formatted_coefficient


def format_counter(counter):
    counter = counter - (counter//100)*100
    formatted_counter = str(counter)
    if len(formatted_counter) == 1:
        formatted_counter = ' '+formatted_counter
    return formatted_counter


def not_found(filename):
    exit("Error: file '"+filename+"' not found!")


def number2binarylist(number):
    temporary = "{0:b}".format(number)
    binary_list = []
    for i in range(len(temporary)-1, -1, -1):
        binary_list.append(int(temporary[i]))
    return binary_list
