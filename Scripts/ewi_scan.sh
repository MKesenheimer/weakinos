#!/bin/bash
#

mkdir ewi_scan
cp -r testrun_clean ./ewi_scan/
cp pwhg_main_nixj merge-pwg-stat runparallel.sh ./ewi_scan/
cd ./ewi_scan/

N=6
for i in `seq 1 $N`; do
  EWI="1d-$i"
  echo "ewi = $EWI"
  ./runparallel.sh -p 4 -g -c -e pwhg_main_nixj -d "run_input_nsusy_n2x1+_ewi$EWI" \
  --fin1 1000023 --fin2 1000024 --slha input_nsusy_1410.4999.slha \
  --ewi $EWI > "log_input_nsusy_n2x1+_ewi$EWI"
  for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
  done
done

for i in `seq 1 $N`; do
  EWI="3d-$i"
  echo "ewi = $EWI"
  ./runparallel.sh -p 4 -g -c -e pwhg_main_nixj -d "run_input_nsusy_n2x1+_ewi$EWI" \
  --fin1 1000023 --fin2 1000024 --slha input_nsusy_1410.4999.slha \
  --ewi $EWI > "log_input_nsusy_n2x1+_ewi$EWI"
  for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
  done
done