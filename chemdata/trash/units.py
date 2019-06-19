"""
"""

from ase import units as ase_units


def kmg_unit(string):
    string = string.lower()
    if 'k' in string:
        return 2**10
    elif 'm' in string:
        return 2**20
    elif 'g' in string:
        return 2**30
    elif 't' in string:
        return 2**40
    elif 'p' in string:
        return 2**50
    else:
        raise NotImplementedError(string)

def unit_to_Ang(unit=None):
    if unit is None:
        return 1
    assert isinstance(unit, str), 'unit should be a string'
    unit = unit.lower()
    if 'ang' in unit:
        return 1
    if 'bohr' in unit or 'au' in unit:
        return ase_units.u['Bohr']
    elif 'nm' in unit:
        return ase_units.u['nm']
    elif 'pm' in unit:
        return ase_units.u['nm']/1000
    else:
        raise NotImplementedError(unit, ' not valid')


