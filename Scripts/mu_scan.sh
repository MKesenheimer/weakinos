#!/bin/bash
# Copyright (C) Matthias Kesenheimer - All Rights Reserved
# Written by Matthias Kesenheimer <m.kesenheimer@gmx.net>, 2017

PWD=$(pwd)
PROC=n1n2
DIR=$PWD/mu_scan_$PROC

mkdir $DIR
cp -r testrun_clean $DIR
cp pwhg_main_ninj merge-pwg-stat merge-data runparallel.sh $DIR
cd $DIR

N=3
for i in `seq -$N 1 $N`; do
for j in `seq -$N 1 $N`; do
  MUR=$(bc -l <<< "scale=3; 2^$i")
  MUR="${MUR}d0"
  MUF=$(bc -l <<< "scale=3; 2^$j")
  MUF="${MUF}d0"
  echo "mur = $MUR, muf = $MUF"
  ./runparallel.sh -p 4 -g -c -e pwhg_main_ninj -d "run_mSUGRA_${PROC}_muR${MUR}_muF${MUF}" \
  --fin1 1000022 --fin2 1000022 --slha input_mSUGRA_1410.4999.slha --merge \
  --mur $MUR --muf $MUF > "log_mSUGRA_${PROC}_muR${MUR}_muF${MUF}"
  for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
  done
done
done

# generate the output file
echo "# muR muF sig err" > mu_scan_results
grep --with-filename "total (btilde" ./log* | \
sed "s/.\\/log_mSUGRA_${PROC}_muR//g" | \
sed "s/: total (btilde+remnants+regulars+osresR) cross section in pb  / /g" | \
sed "s/  +-    / /g" | \
sed "s/_muF/  /g">> mu_scan_results