"""
independent chemical symbols
"""


__version__ = '1.0.0'

import os
BASEDIR = os.path.dirname(os.path.abspath(__file__))
atom_file = os.path.join(BASEDIR, 'atom.json')

chemical_symbols = [
    # 0
    'X',
    # 1
    'H', 'He',
    # 2
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    # 3
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    # 4
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    # 5
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    # 6
    'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
    'Po', 'At', 'Rn',
    # 7
    'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
    'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc',
    'Lv', 'Ts', 'Og']

def version():
    return __version__


global atom_dict
atom_dict = None
def read_atom():
    global atom_dict
    if atom_dict is None:
        assert os.path.exists(atom_file), 'atom file does not exist'
        with open(atom_file) as fd:
            import json
            atom_dict = json.load(fd)


def get_element_property(element, _property='covalentradius'):
    read_atom()
    return atom_dict[element][_property]


def get_element_index(element):
    return chemical_symbols.index(element)


def get_element_covalent(element):
    read_atom()
    return atom_dict[element]['covalentradius']

def get_element_mass(element):
    read_atom()
    return atom_dict[element]['mass']


def get_element_electronic_configuration(element):
    read_atom()
    return atom_dict[element]['Electron_Configuration']

def get_all_properties():
    read_atom()
    return list(atom_dict['X'].keys())

