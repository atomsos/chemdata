#!/usr/bin/env python3

import os, sys
import argparse

from chemdata import chemical_symbols
from itertools import combinations



def main(min_length, max_length, connector):
    for length in range(min_length, max_length+1):
        for syms in combinations(chemical_symbols[1:], length):
            print(connector.join(syms))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("min_length_max_length", nargs=2, metavar='int')
    parser.add_argument("connector", nargs='?', default='-')
    args = parser.parse_args()
    # print(args)
    min_length, max_length = [int(x) for x in args.min_length_max_length]
    main(min_length, max_length, args.connector)

