#!/bin/bash


CURDIR="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")" 

export PYTHONPATH=$CURDIR:$PYTHONPATH
export PATH=$CURDIR/bin:$PATH
