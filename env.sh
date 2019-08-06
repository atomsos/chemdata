#!/bin/bash


[ "$(uname)" == "Linux" ] && CURDIR="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")" 
[ "$(uname)" == "Darwin" ] && CURDIR="$(dirname -- "$(readlink -- -f "${BASH_SOURCE[0]}")")"
[ "$(uname)" == "Darwin" ] && [ "$CURDIR" == "." ] && CURDIR="$(dirname -- "${PWD}/${BASH_SOURCE[0]}")"

CURDIR=$(cd $CURDIR; pwd)
echo $CURDIR


export PYTHONPATH=$CURDIR:$PYTHONPATH
export PATH=$CURDIR/bin:$PATH
