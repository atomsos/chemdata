"""
test chemdata
"""


import chemdata

    
def test():
    print(chemdata.version())
    properties = chemdata.get_all_properties()
    for sym in chemdata.chemical_symbols[:100]:
        for p in properties:
            chemdata.get_element_property(sym, p)
    for func_name in chemdata.__dir__():
        if func_name.startswith('get_element_'):
            data = getattr(chemdata, func_name)('H')
            print(func_name, data, type(data))
        elif func_name.startswith('get_spectra_'):
            data = getattr(chemdata, func_name)('V', 10)
            print(func_name, data, type(data))


if __name__ == '__main__':
    test()
