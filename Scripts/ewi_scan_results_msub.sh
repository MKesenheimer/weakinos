#!/bin/bash
# Copyright (C) Matthias Kesenheimer - All Rights Reserved
# Written by Matthias Kesenheimer <m.kesenheimer@gmx.net>, 2017

PWD=$(pwd)
PROC=n1n1
DIR=$PWD/ewi_scan_$PROC

cd $DIR

rm $DIR/ewi_scan_results_${PROC}
echo "# ewi sig err" > $DIR/ewi_scan_results_${PROC}

N=6
for i in `seq 1 $N`; do
  EWI="1d-$i"
  RUNDIR=${DIR}/run_input_orig_${PROC}_ewi${EWI}
  cd $RUNDIR
  rm pwg-st2-combined-stat.dat
  ../merge-pwg-stat $(ls ./pwg-st2-*-stat.dat) > pwg-st2-combined-stat.dat
  grep "total (btilde" pwg-st2-combined-stat.dat | \
  sed "s/ total (btilde+remnants+regulars+osresR) cross section in pb    /$EWI /g" | \
  sed "s/  +-   / /g" | tee -a $DIR/ewi_scan_results_${PROC}
done

for i in `seq 1 $N`; do
  EWI="3d-$i"
  RUNDIR=${DIR}/run_input_orig_${PROC}_ewi${EWI}
  cd $RUNDIR
  rm pwg-st2-combined-stat.dat
  ../merge-pwg-stat $(ls ./pwg-st2-*-stat.dat) > pwg-st2-combined-stat.dat
  grep "total (btilde" pwg-st2-combined-stat.dat | \
  sed "s/ total (btilde+remnants+regulars+osresR) cross section in pb    /$EWI /g" | \
  sed "s/  +-   / /g" | tee -a $DIR/ewi_scan_results_${PROC}
done