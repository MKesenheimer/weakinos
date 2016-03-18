#!/bin/bash
#

PWD=$(pwd)
PROC=n1n1
DIR=$PWD/mu_scan_$PROC

cd $DIR

# generate the output file
echo "# muR muF sig err" > mu_scan_results_${PROC}
grep --with-filename "total (btilde" ./log* | \
sed "s/.\\/log_input_mSUGRA_${PROC}_muR//g" | \
sed "s/: total (btilde+remnants+regulars+osresR) cross section in pb  / /g" | \
sed "s/  +-   / /g" | \
sed "s/_muF/  /g">> mu_scan_results_${PROC}
