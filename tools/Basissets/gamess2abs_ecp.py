import sys
if sys.version_info[0] < 3:
    sys.exit('This script requires Python 3')

import os
import re

elem_symbs = (
    'X',
    'H', 'He',
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga',
        'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
        'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
        'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt',
        'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
    'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
        'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
        'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og')
elem_names = (
    '',
    'Hydrogen', 'Helium',
    'Lithium', 'Beryllium', 'Boron', 'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine',
        'Neon',
    'Sodium', 'Magnesium', 'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur',
        'Chlorine', 'Argon',
    'Potassium', 'Calcium', 'Scandium', 'Titanium', 'Vanadium', 'Chromium',
        'Manganese', 'Iron', 'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium',
        'Germanium', 'Arsenic', 'Selenium', 'Bromine', 'Krypton',
    'Rubidium', 'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum',
        'Technetium', 'Ruthenium', 'Rhodium', 'Palladium', 'Silver', 'Cadmium',
        'Indium', 'Tin', 'Antimony', 'Tellurium', 'Iodine', 'Xenon',
    'Caesium', 'Barium', 'Lanthanum', 'Cerium', 'Praseodymium', 'Neodymium',
        'Promethium', 'Samarium', 'Europium', 'Gadolinium', 'Terbium',
        'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium', 'Lutetium',
        'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium', 'Osmium', 'Iridium',
        'Platinum', 'Gold', 'Mercury', 'Thallium', 'Lead', 'Bismuth',
        'Polonium', 'Astatine', 'Radon',
    'Francium', 'Radium', 'Actinium', 'Thorium', 'Protactinium', 'Uranium',
        'Neptunium', 'Plutonium', 'Americium', 'Curium', 'Berkelium',
        'Californium', 'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium',
        'Lawrencium', 'Rutherfordium', 'Dubnium', 'Seaborgium', 'Bohrium',
        'Hassium', 'Meitnerium', 'Darmstadtium', 'Roentgenium', 'Copernicium',
        'Nihonium', 'Flerovium', 'Moscovium', 'Livermorium', 'Tennessine',
        'Oganesson')
AO_types = 'SPDFGHIKLMNOQRTUVWXYZ'
is_arr = lambda obj: hasattr(obj, '__len__') and hasattr(obj, '__getitem__') \
    and not isinstance(obj, str)

def read(filename):
    with open(filename, 'r') as file:
        string = file.read()
    return string


def write(filename, string):
    with open(filename, 'w+') as file:
        file.write(string)


def search(pattern, string, flags = 0):
    res = re.search(pattern, string, flags)
    if res:
        res = res.groups()
        if len(res) == 1:
            res = res[0]
        return res


def findall(pattern, string, flags = 0):
    return re.findall(pattern, string, flags)


def sub(patterns, string, flags = 0):
    if not is_arr(patterns[0]):
        patterns = (patterns,)
    for old_pattern, new_pattern in patterns:
        string = re.sub(old_pattern, new_pattern, string, flags = flags)
    return string


def str_between(string, after, before, idx = 0):
    start = 0
    while idx > -1:
        idx -= 1
        new_start = string.find(after, start)
        if new_start == -1:
            break
        else:
            start = new_start + len(after)
    end = string.find(before, start)
    if end == -1:
        end = len(string)
    return string[start : end]

if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("gamess2abs_ecp.py converts a basis set file with ECPs " \
            + "in GAMESS US format (from Basis Set Exchange) to basis set " \
            + "files in Amolqc format ('.abs' and '.ecp')\n" \
            + "usage:\n" \
            + "gamess2abs_ecp.py basis_file\n")
    filename = sys.argv[1]
    name = os.path.splitext(filename)[0]
        
    text = read(filename)
    comment_text = search('((?:!.+\n)+)', text)
    bas_text = str_between(text, '$DATA', '$END')
    ecp_text = str_between(text, '$ECP', '$END')

    separator = '****\n'

    text = bas_text
    for name in findall(r'([A-Z]{3,})', text):
        text = text.replace('\n' + name, separator \
            + elem_symbs[elem_names.index(name[0] + name[1 :].lower())], 1)
    text = sub((r'\n\d+', '\n'), text)
    text = comment_text.replace('!', '#') \
        + sub((r'\b([' + AO_types + ']) +(\d+)', r'\g<1> GTO \g<2>'), \
        text)[len(separator) :] + separator
    write(filename + '.abs', text)
    
    text = ecp_text
    for substr, symb in findall(r'((\w{1,2})-ECP GEN)', text):
        if len(symb) == 2:
            symb = symb[0] + symb[1].lower()
        text = text.replace(substr, ' ' + separator + '####' \
            + elem_names[elem_symbs.index(symb)] + '\n' \
            + str(elem_symbs.index(symb)), 1)
    text = comment_text.replace('!', '#') \
        + sub((r' *-{5}.+? potential -{5}', ''), text)[len(separator) + 1 :] \
        + ' ' + separator
    write(filename + '.ecp', text)
