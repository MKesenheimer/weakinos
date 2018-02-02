#!/bin/bash
# Copyright (C) Matthias Kesenheimer - All Rights Reserved
# Written by Matthias Kesenheimer <m.kesenheimer@gmx.net>, 2017

PWD=$(pwd)
PROC=n1n1
DIR=$PWD/ewi_scan_$PROC

mkdir $DIR
cp -r testrun_clean $DIR
cp pwhg_main_ninj merge-pwg-stat merge-data runparallel.sh $DIR
cd $DIR

N=6
for i in `seq 1 $N`; do
  EWI="1d-$i"
  echo "submitting ewi = $EWI"
  ./runparallel.sh -p 10 -g -c -e pwhg_main_ninj -d "run_input_orig_${PROC}_ewi${EWI}" \
  --fin1 1000022 --fin2 1000022 --slha input_orig.slha \
  --ncall1 200000 --ncall2 300000 --ncall1osres 200000 --ncall2osres 300000 \
  --ewi $EWI \
  --usemsub --offset 0 > "log_input_orig_${PROC}_ewi$EWI"
done

for i in `seq 1 $N`; do
  EWI="3d-$i"
  echo "submitting ewi = $EWI"
  ./runparallel.sh -p 10 -g -c -e pwhg_main_ninj -d "run_input_orig_${PROC}_ewi${EWI}" \
  --fin1 1000022 --fin2 1000022 --slha input_orig.slha \
  --ncall1 200000 --ncall2 300000 --ncall1osres 2000000 --ncall2osres 300000 \
  --ewi $EWI \
  --usemsub --offset 0 > "log_input_orig_${PROC}_ewi$EWI"
done

# generate the output file
#echo "# ewi sig err" > ewi_scan_results
#grep --with-filename "total (btilde" ./log* | \
#sed "s/.\\/log_input_nsusy_${PROC}_ewi//g" | \
#sed "s/: total (btilde+remnants+regulars+osresR) cross section in pb  / /g" | \
#sed "s/  +-    / /g" >> ewi_scan_results