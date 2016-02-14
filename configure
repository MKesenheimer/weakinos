#!/bin/sh
# calls the LoopTools and SLHAlib configure scripts

WORKINGDIR=$(pwd)
LT=$WORKINGDIR/Tools/LoopTools-2.12
SLHA=$WORKINGDIR/Tools/SLHALib-2.2
PYTHIA=$WORKINGDIR/Tools/pythia8215

#cd $PYTHIA && ./configure --prefix-lib=$WORKINGDIR/Tools --cxx-common='-O2 -fomit-frame-pointer -ffast-math -Wall -m64 -stdlib=libstdc++ -mmacosx-version-min=10.6 -Qunused-arguments -g -ansi -pedantic -W -Wall -Wshadow -fPIC'
cd $PYTHIA && ./configure --prefix-lib=$WORKINGDIR/Tools
cd $LT && ./configure
cd $SLHA && ./configure