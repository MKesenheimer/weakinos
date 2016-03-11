#!/bin/bash
#

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
  echo "ewi = $EWI"
  ./runparallel.sh -p 4 -g -c -e pwhg_main_ninj -d "run_input_nsusy_${PROC}_ewi${EWI}" \
  --fin1 1000022 --fin2 1000022 --slha input_nsusy_1307.0782.slha --merge \
  --ewi $EWI > "log_input_nsusy_${PROC}_ewi$EWI"
  for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
  done
done

for i in `seq 1 $N`; do
  EWI="3d-$i"
  echo "ewi = $EWI"
  ./runparallel.sh -p 4 -g -c -e pwhg_main_ninj -d "run_input_nsusy_${PROC}_ewi$EWI" \
  --fin1 1000022 --fin2 1000022 --slha input_nsusy_1307.0782.slha --merge \
  --ewi $EWI > "log_input_nsusy_${PROC}_ewi$EWI"
  for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
  done
done

# generate the output file
echo "# ewi sig err" > ewi_scan_results
grep --with-filename "total (btilde" ./log* | \
sed "s/.\\/log_input_nsusy_${PROC}_ewi//g" | \
sed "s/: total (btilde+remnants+regulars+osresR) cross section in pb  / /g" | \
sed "s/  +-    / /g" >> ewi_scan_results