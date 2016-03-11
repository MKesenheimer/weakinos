#!/bin/bash
#

PWD=$(pwd)
PROC=n1n2
DIR=$PWD/mu_scan_$PROC

#mkdir $DIR
#cp -r testrun_clean $DIR
#cp pwhg_main_ninj merge-pwg-stat merge-data runparallel.sh $DIR
cd $DIR

N=3
for i in `seq -$N 1 $N`; do
for j in `seq $N -1 -$N`; do
  MUR=$(bc -l <<< "scale=1; 2^$i")
  MUR="${MUR}d0"
  MUF=$(bc -l <<< "scale=1; 2^$j")
  MUF="${MUF}d0"
  MUR2=$(bc -l <<< "scale=3; 2^$i")
  MUR2="${MUR2}d0"
  MUF2=$(bc -l <<< "scale=3; 2^$j")
  MUF2="${MUF2}d0"
  echo "mur = $MUR, muf = $MUF"
  RUNDIR=run_input_nsusy_${PROC}_mur${MUR}_muF${MUF}
  #$DIR/merge-pwg-stat $(ls $RUNDIR/pwg-st2-*-stat.dat) > $DIR/log_input_nsusy_${PROC}_muR${MUR}_muF${MUF}
  cat $RUNDIR/pwg-st2-combined-stat.dat > $DIR/log_input_nsusy_${PROC}_muR${MUR2}_muF${MUF2}
done
done

# generate the output file
echo "# muR muF sig err" > mu_scan_results
grep --with-filename "total (btilde" ./log* | \
sed "s/.\\/log_input_nsusy_${PROC}_muR//g" | \
sed "s/: total (btilde+remnants+regulars+osresR) cross section in pb  / /g" | \
sed "s/  +-   / /g" | \
sed "s/_muF/  /g">> mu_scan_results
