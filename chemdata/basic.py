"""
independent chemical symbols
"""




__all__ = [
    'get_element',
    'get_element_number',
    'get_element_property',
    'get_element_covalent',
    'get_element_mass',
    'get_element_electronic_configuration',
    'get_element_pair_covalent',
    'get_all_properties',
    'chemical_symbols',
]




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

def get_element(element):
    if isinstance(element, str):
        return element
    else:
        try:
            element = int(element)
        except Exception as e:
            print(element, 'is not a symbol or int')
            raise e
        element = chemical_symbols[element]
    return element


def get_bash_elements():
    print(' '.join(chemical_symbols[1:]))

global ATOM_DICT
ATOM_DICT = None
def READ_ATOM():
    global ATOM_DICT
    if ATOM_DICT is None:
        assert os.path.exists(atom_file), 'atom file does not exist'
        with open(atom_file) as fd:
            import json
            ATOM_DICT = json.load(fd)


def get_element_number(element):
    element = get_element(element)
    return chemical_symbols.index(element)


def get_element_property(element, _property='mass'):
    READ_ATOM()
    element = get_element(element)
    return ATOM_DICT[element][_property]


def get_all_properties():
    READ_ATOM()
    return list(ATOM_DICT['X'].keys())


def get_element_covalent(element):
    return get_element_property(element, 'covalentradius')


def get_element_mass(element):
    return get_element_property(element, 'mass')

def get_element_electronic_configuration(element):
    return get_element_property(element, 'Electron_Configuration')


def get_valence_electron(element):
    READ_ATOM()
    config = get_element_electronic_configuration(element)
    full_shell = {
        's' : 2,
        'p' : 6,
        'd' : 10,
        'f' : 14,
        'g' : 16,
    }
    is_full = True
    n_valence = 0
    for cfg in config.split('_'):
        if cfg.startswith('['):
            continue
        n, l = cfg[0], cfg[1]
        nele = int(cfg[2:])
        if is_full and nele < full_shell[l]:
            is_full = False
        if not is_full:
            n_valence += nele
    if n_valence == 0:
        n_valence = nele
    return n_valence

def get_element_pair_covalent(element1, element2=None):
    return get_element_covalent(element1) + get_element_covalent(element2 or element1)
