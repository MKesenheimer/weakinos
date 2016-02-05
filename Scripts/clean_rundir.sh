#!/bin/bash

WORKINGDIR=$(pwd)
RUNDIR=$WORKINGDIR/$1

find $RUNDIR ! \( -name '*.slha' -o -name '*.input' -o -name 'pwgseeds.dat' \) -type f -exec rm -f {} +
