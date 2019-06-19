"""
NIST data
https://physics.nist.gov/PhysRefData/ASD/lines_form.html
"""


__all__ = [
    'get_spectra_name',
    'get_spectra_ground_shells',
    'get_spectra_ground_level',
    'get_spectra_ionized_level',
    'get_spectra_ionization_energy',
    'get_spectra_uncertainty',
    'get_spectra_references',
]

import os
import re
import json
from . import basic

BASEDIR = os.path.dirname(os.path.abspath(__file__))
NIST_SPECTRUM_FILE = os.path.join(BASEDIR, 'nist_spectrum.json')


global NIST_SPECTURM
NIST_SPECTURM = None


def LOAD_NIST_SPECTURM():
    global NIST_SPECTURM
    if NIST_SPECTURM is None:
        with open(NIST_SPECTRUM_FILE) as fd:
            NIST_SPECTURM = json.load(fd)


def PARSE_DATA(data):
    float_pattern = r'^[\(\[]?(\d+\s*\d*\.?\d*)[\)\]]?$'
    if re.match(float_pattern, data):
        data = float(re.match(float_pattern, data)[1])
    return data


def GET_DATA(element, charge, key):
    global NIST_SPECTURM
    LOAD_NIST_SPECTURM()
    # import pdb; pdb.set_trace()
    element_number = str(basic.get_element_number(element))
    if not element_number in NIST_SPECTURM:
        raise ValueError('element: ', element, 'not in NIST_SPECTURM')
    NIST_SPECTURM_ELEMENT = NIST_SPECTURM[element_number]
    charge = str(int(charge))
    if not charge in NIST_SPECTURM_ELEMENT:
        raise ValueError('charge: ', charge, 'not in NIST_SPECTURM for', element)
    DATA = NIST_SPECTURM_ELEMENT[charge]
    for dkey, dval in DATA.items():
        if dkey.replace(' ', '').lower().startswith(key.lower().replace(' ', '')):
            if dval:
                return PARSE_DATA(dval)
            return None
    raise ValueError(key, 'not exist for element', element)


def get_spectra_name(element, charge):
    return GET_DATA(element, charge, 'Sp. Name')

def get_spectra_ground_shells(element, charge):
    return GET_DATA(element, charge, 'Ground Shells')

def get_spectra_ground_level(element, charge):
    return GET_DATA(element, charge, 'Ground Level')

def get_spectra_ionized_level(element, charge):
    return GET_DATA(element, charge, 'Ionized Level')

def get_spectra_ionization_energy(element, charge):
    return GET_DATA(element, charge, 'Ionization Energy')

def get_spectra_uncertainty(element, charge):
    return GET_DATA(element, charge, 'Uncertainty')

def get_spectra_references(element, charge):
    return GET_DATA(element, charge, 'References')

