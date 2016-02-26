#!/bin/bash
#

mkdir mu_scan
cp -r testrun_clean ./mu_scan/
cp pwhg_main_nixj merge-pwg-stat runparallel.sh ./mu_scan/
cd ./mu_scan/

N=1
for i in `seq -$N 1 $N`; do
for j in `seq -$N 1 $N`; do
  MUR=$(bc -l <<< "scale=1; 2^$i")
  MUR="${MUR}d0"
  MUF=$(bc -l <<< "scale=1; 2^$j")
  MUF="${MUF}d0"
  echo "mur = $MUR, muf = $MUF"
  ./runparallel.sh -p 4 -g -c -e pwhg_main_nixj -d "run_input_nsusy_wevents_n2x1+_mur${MUR}_muF${MUF}" \
  --fin1 1000023 --fin2 1000024 --slha input_nsusy_1410.4999.slha --genevents \
  --mur $MUR --muf $MUF > "log_input_nsusy_wevents_n2x1+_mur${MUR}_muF${MUF}"
  for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
  done
done
done