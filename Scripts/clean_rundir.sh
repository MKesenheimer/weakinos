#!/bin/bash
# Copyright (C) Matthias Kesenheimer - All Rights Reserved
# Written by Matthias Kesenheimer <m.kesenheimer@gmx.net>, 2017

WORKINGDIR=$(pwd)
RUNDIR=$WORKINGDIR/$1

if [[ $1 == "" ]]; then
    echo "usage: ./clean_rundir.sh <directory>"
    exit -1
fi

find $RUNDIR ! \( -name '*.slha' -o -name '*.input' -o -name 'pwgseeds.dat' \) -type f -exec rm -f {} +
