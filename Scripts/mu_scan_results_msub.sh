#!/bin/bash
#

PWD=$(pwd)
PROC=n1n2
DIR=$PWD/mu_scan_$PROC

cd $DIR

# generate the output file
echo "Combined results for stage 2:"
echo "muR muF sig err"
echo "# muR muF sig err" > $DIR/mu_scan_results_${PROC}

N=3
for i in `seq -$N 1 $N`; do
for j in `seq -$N 1 $N`; do
  MUR=$(bc -l <<< "scale=3; 2^$i")
  MUR="${MUR}d0"
  MUF=$(bc -l <<< "scale=3; 2^$j")
  MUF="${MUF}d0"
  RUNDIR=${DIR}/run_mSUGRA_${PROC}_mur${MUR}_muF${MUF}
  cd $RUNDIR
  rm pwg-st2-combined-stat.dat
  ../merge-pwg-stat $(ls ./pwg-st2-*-stat.dat) > pwg-st2-combined-stat.dat
  grep "total (btilde" pwg-st2-combined-stat.dat | \
  sed "s/ total (btilde+remnants+regulars+osresR) cross section in pb    /$MUR $MUF /g" | \
  sed "s/  +-   / /g" | tee -a $DIR/mu_scan_results_${PROC}
done
done