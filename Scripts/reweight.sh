#!/bin/bash
# Copyright (C) Matthias Kesenheimer - All Rights Reserved
# Written by Matthias Kesenheimer <m.kesenheimer@gmx.net>, 2017
#
# Examples:
# ./reweight.sh -e pwhg_main_nixj -d testrun1 -r 1 --muf 0.5 --merge
# -> reweight the events in testrun1 with index = 1 and muf = 0.5. At the
#    end merge the event files.

# Functions
function info {
   echo "Reweights the events in parallel mode"
}

function usage {
cat <<EOM
Usage: $(basename $0) [OPTIONS]
Mandatory arguments:
  -d, --directory <name>   specify the subdirectory relative to working
                           directory where to execute the program and where 
                           to look for powheg.input
  -e, --executable <name>  the relative path of the executable
  -r, --lhrwgtid           the id of the reweighting

Optional arguments:
  -p, --parallel <n>       number of parallel jobs to submit
  -o, --offset <n>         supply an offset for extracting the seeds from seeds.dat
  --mur <n>                set renscfact
  --muf <n>                set facscfact
  --pdf <n>                LHA number of pdf
  --merge                  merge the event files and delete the old ones
EOM
   exit 0
}

function overwrite_powheg_var {
  sed -i -e "s/$1 /#$1 /g" powheg.input
  echo "$1 $2" >> powheg.input
}

function comment_powheg_var {
  sed -i -e "s/$1 /#$1 /g" powheg.input
}

# make sure grep uses no other arguments (for example via aliases)
function read_powheg_var {
  grep "$1" powheg.input | sed 's/[^0-9]//g' | head -1
}

# if there are no arguments print help
if [ $# -lt 1 ] 
  then usage
fi

#default values
JOBS=4
NICENESS=10

# go through the options
while [[ $# -gt 0 ]]; do
KEY="$1"
case $KEY in
    -h|--help)
        usage
        shift # past argument
        ;;
    -i|--info)
        info
        shift
        ;;
    -d|--directory)
        RUNDIR="$2"
        shift
        shift
        ;;
    -e|--executable)
        EXE="$2"
        shift
        shift
        ;;
    -r|--lhrwgtid)
        ID="$2"
        shift
        shift
        ;;
    -p|--parallel)
        JOBS="$2"
        shift
        shift
        ;;
    -o|--offset)
        NSEEDOFFSET="$2"
        shift
        shift
        ;;
    --muf)
        FACSCFACT="$2"
        shift
        shift
        ;;
    --mur)
        RENSCFACT="$2"
        shift
        shift
        ;;
    --pdf)
        LOPDF="$2"
        shift
        shift
        ;;
    --merge)
        MERGE=true
        shift
        ;;
    *)
        usage    # unknown option
        ;;
esac
done

# check if RUNDIR is set
if [ "$RUNDIR" == "" ]; then
   echo "Error: no directory specified."
   usage
fi

# check if EXE is set
if [ "$EXE" == "" ]
then
   echo "Error: no executable specified."
   usage
fi

# check if ID is set
if [ "$ID" == "" ]
then
   echo "Error: no reweighting id specified."
   usage
fi

# directories
WORKINGDIR=${PWD}
RUNDIR=$WORKINGDIR/$RUNDIR
EXEPATH=$WORKINGDIR/$EXE

# change into the directory where to execute pwhg_main
cd $RUNDIR

echo ""
echo "Reweighting the events"
echo "" >> powheg.input
echo "# Reweighting the events" >> powheg.input

ID=$(printf "%03d\n" $ID)
overwrite_powheg_var "lhrwgt_id" \'"$ID"\'
overwrite_powheg_var "compute_rwgt" 1
comment_powheg_var "testplots"

# sed magic
if [ "$RENSCFACT" != "" ]; then
   overwrite_powheg_var "renscfact" "$RENSCFACT"
fi

if [ "$FACSCFACT" != "" ]; then
   overwrite_powheg_var "facscfact" "$FACSCFACT"
fi

if [ "$PDF" != "" ]; then
   overwrite_powheg_var "lhans1" $PDF
   overwrite_powheg_var "lhans2" $PDF
fi

# if reweighting more than once -> rename reweighted files to pwgevents-* and start again
if ls pwgevents-rwgt-* 1> /dev/null 2>&1; then
  # files do exist
  find $RUNDIR -name "pwgevents-rwgt-*" -exec bash -c 'mv $0 ${0/pwgevents-rwgt/pwgevents}' {} \;
fi

# start the POWHEG-main executable
echo "  starting $JOBS job(s)..."
for i in `seq 1 $JOBS`; do
   NSEED=$((i+NSEEDOFFSET))
   echo "  job $i with nseed $NSEED"
   INDEX=$(printf "%04d\n" $NSEED)
   nohup nice -n $NICENESS $EXEPATH < <(printf "%s\n" "$NSEED" "pwgevents-${INDEX}.lhe") > $RUNDIR/powheg_reweight_$i.output 2>&1 &
done

for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
done

if [ "$MERGE" = true ]; then
  # merge the event files
  cat $RUNDIR/pwgevents-rwgt-*.lhe | grep -v "/LesHouchesEvents" > $RUNDIR/pwgevents-rwgt.lhe
  echo "</LesHouchesEvents>" >> $RUNDIR/pwgevents-rwgt.lhe
fi  