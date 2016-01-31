#!/bin/bash
#
# examples:
# in ./weakinos:
# ./runparallel.sh -d ./neuIchaJ/testrun_1 -e ./neuIchaJ/pwhg_main_nixj -p 4
# -> run pwhg_main_nixj in testrun_1 on 4 cores
#
# in ./weakinos/neuIchaJ:
# ../runparallel.sh -c -d ./testrun_1 -e ./pwhg_main_nixj -p 4 --itmx1 4 \
# --itmx2 4 --itmx1osres 6 --itmx2osres 8 --ncall1 2000 --ncall2 2000 \
# --ncall1osres 20000 --ncall2osres 20000 -o 0
# -> run pwhg_main_nixj in testrun_1 on 4 cores and overwrite some powheg
#    parameters in powheg.input
#
# automatically runs POWHEG in parallel mode with the following steps:
#Step 1a:
# - fakevirtuals 1
# - xgriditeration 1
# - parallelstage 1
#
#Step 1b:
# - fakevirtuals 1
# - xgriditeration 2
# - parallelstage 1
#
#Step 2: 
# - fakevirtuals 0
# - xgriditeration 1
# - parallelstage 2
# - numevts 0
# - nubound 0
#
#Additionally:
#Step 3:
# - nubound 100000
# - parallelstage 3
#
#Step 4:
# - numevts 100000
# - parallelstage 4

# Functions
function info {
   echo "Automatically runs POWHEG in parallel mode"
}

function usage {
cat <<EOM
Usage: $(basename $0) [OPTIONS]
Mandatory arguments:
  -d, --directory <name>   specify the subdirectory relative to working
                           directory where to execute the program and where 
                           to look for powheg.input
  -e, --executable <name>  the relative path of the executable
  -p, --parallel <n>       number of parallel jobs to submit
  
Optional arguments:
  -h, --help               print this help message
  -i, --info               show informations
  --ncall1 <n>             overwrite the parameter ncall1 in powheg.input
  --ncall2 <n>             overwrite the parameter ncall2 in powheg.input
  --ncall1osres <n>        overwrite the parameter ncall1osres in powheg.input
  --ncall2osres <n>        overwrite the parameter ncall1osres in powheg.input
  --itmx1 <n>              overwrite the parameter itmx1 in powheg.input
  --itmx2 <n>              overwrite the parameter itmx2 in powheg.input
  --itmx1osres <n>         overwrite the parameter itmx1osres in powheg.input
  --itmx2osres <n>         overwrite the parameter itmx2osres in powheg.input
  -c, --clean              clean directory before running POWHEG
  -o, --offset <n>         supply an offset for extracting the seeds from seeds.dat
  --genevents <n1> <n2>    generate the upper bound (n1) and events (n2)
  --usemsub                use the submitting system msub
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

# if there are no arguments print help
if [ $# -lt 1 ] 
  then usage
fi

#default values
JOBS=1
CLEAN=false
NSEEDOFFSET=0
NICENESS=10
USEMSUB=false

# go through the options
while [[ $# -gt 0 ]]
do
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
    -c|--clean)
        CLEAN=true
        shift
        ;;
    --ncall1)
        NCALL1="$2"
        shift
        shift
        ;;
    --ncall2)
        NCALL2="$2"
        shift
        shift
        ;;
    --ncall1osres)
        NCALL1OSRES="$2"
        shift
        shift
        ;;
    --ncall2osres)
        NCALL2OSRES="$2"
        shift
        shift
        ;;
    --itmx1)
        ITMX1="$2"
        shift
        shift
        ;;
    --itmx2)
        ITMX2="$2"
        shift
        shift
        ;;
    --itmx1osres)
        ITMX1OSRES="$2"
        shift
        shift
        ;;
    --itmx2osres)
        ITMX2OSRES="$2"
        shift
        shift
        ;;
    --genevents)
        NEVENTS="$2"
        NUBOUND="$3"
        shift
        shift
        shift
        ;;
    -usemsub)
        USEMSUB=true
        shift
        ;;
    --default)
        DEFAULT=YES
        ;;
    *)
        usage    # unknown option
        ;;
esac
done

# check if RUNDIR is set
if [ "$RUNDIR" == "" ]
then
   echo "Error: no directory specified."
   usage
fi

# check if EXE is set
if [ "$EXE" == "" ]
then
   echo "Error: no executable specified."
   usage
fi

# directories
WORKINGDIR=${PWD}
RUNDIR=$WORKINGDIR/$RUNDIR
EXEPATH=$WORKINGDIR/$EXE

# change into the directory where to execute pwhg_main
cd $RUNDIR

# clean up the directory
if [ "$CLEAN" = true ]
then
   find $RUNDIR ! \( -name '*.slha' -o -name '*.input' -o -name 'pwgseeds.dat' \) -type f -exec rm -f {} +
   cp powheg_clean.input powheg.input
fi
#exit 0

# append to powheg.input
echo "" >> powheg.input
echo "" >> powheg.input
echo "# Modified by ../runparallel:" >> powheg.input

#default parameters in powheg.input
overwrite_powheg_var "use-old-grid" 1
overwrite_powheg_var "use-old-ubound" 1
comment_powheg_var "iseed"
overwrite_powheg_var "manyseeds" 1
overwrite_powheg_var "testplots" 1

# sed magic
if [ "$NCALL1" != "" ]
then
   overwrite_powheg_var "ncall1" $NCALL1
fi  

if [ "$NCALL2" != "" ]
then
   overwrite_powheg_var "ncall2" $NCALL2
fi

if [ "$NCALL1OSRES" != "" ]
then
   overwrite_powheg_var "ncall1osres" $NCALL1OSRES
fi  

if [ "$NCALL2OSRES" != "" ]
then
   overwrite_powheg_var "ncall2osres" $NCALL2OSRES
fi

if [ "$ITMX1" != "" ]
then
   overwrite_powheg_var "itmx1" $ITMX1
fi  

if [ "$ITMX2" != "" ]
then
   overwrite_powheg_var "itmx2" $ITMX2
fi

if [ "$ITMX1OSRES" != "" ]
then
   overwrite_powheg_var "itmx1osres" $ITMX1OSRES
fi  

if [ "$ITMX2OSRES" != "" ]
then
   overwrite_powheg_var "itmx2osres" $ITMX2OSRES
fi

# start the POWHEG-main executable

# STEP 1a
echo ""
echo "Step 1a: Generating Grids"
echo "" >> powheg.input
echo "#Step 1a: Generating Grids" >> powheg.input

overwrite_powheg_var "fakevirtuals" 1
overwrite_powheg_var "parallelstage" 1
overwrite_powheg_var "xgriditeration" 1

echo "  starting $JOBS job(s)..."
for i in `seq 1 $JOBS`; do
   NSEED=$((i+NSEEDOFFSET))
   echo "  job $i with nseed $NSEED"
   nohup nice -n $NICENESS $EXEPATH <<< $NSEED > $RUNDIR/powheg_step1a_$i.output 2>&1 &
done

for job in `jobs -p`
do
    wait $job
    echo "  job with pid=$job finished"
done

# STEP 1b
echo ""
echo "Step 1b: Generating Grids"
echo "" >> powheg.input
echo "#Step 1b: Generating Grids" >> powheg.input

overwrite_powheg_var "fakevirtuals" 1
overwrite_powheg_var "parallelstage" 1
overwrite_powheg_var "xgriditeration" 2

echo "  starting $JOBS job(s)..."
for i in `seq 1 $JOBS`; do
   NSEED=$((i+NSEEDOFFSET))
   echo "  job $i with nseed $NSEED"
   nohup nice -n $NICENESS $EXEPATH <<< $NSEED > $RUNDIR/powheg_step1b_$i.output 2>&1 &
done

for job in `jobs -p`
do
    wait $job
    echo "  job with pid=$job finished"
done


# STEP 2
echo ""
echo "Step 2: NLO run"
echo "" >> powheg.input
echo "#Step 2: NLO run" >> powheg.input

overwrite_powheg_var "fakevirtuals" 0
# TODO: xgriditeration = 1 necessary?
overwrite_powheg_var "xgriditeration" 1
overwrite_powheg_var "parallelstage" 2
overwrite_powheg_var "numevts" 0
overwrite_powheg_var "nubound" 0

echo "  starting $JOBS job(s)..."
for i in `seq 1 $JOBS`; do
   NSEED=$((i+NSEEDOFFSET))
   echo "  job $i with nseed $NSEED"
   nohup nice -n $NICENESS $EXEPATH <<< $NSEED > $RUNDIR/powheg_step2_$i.output 2>&1 &
done

for job in `jobs -p`
do
    wait $job
    echo "  job with pid=$job finished"
done

# combined results for stage 2
echo ""
echo "Combined results:"
$WORKINGDIR/merge-pwg-stat $(ls ./pwg-st2-*-stat.dat)

# if the user wants to generate events
if [ "$NEVENTS" != "" ]
then
  # STEP 3
  echo ""
  echo "Step 3: Upper bound"
  echo "" >> powheg.input
  echo "#Step 3: Upper bound" >> powheg.input

  overwrite_powheg_var "nubound" $NUBOUND
  overwrite_powheg_var "parallelstage" 3

  echo "  starting $JOBS job(s)..."
  for i in `seq 1 $JOBS`; do
     NSEED=$((i+NSEEDOFFSET))
     echo "  job $i with nseed $NSEED"
     nohup nice -n $NICENESS $EXEPATH <<< $NSEED > $RUNDIR/powheg_step3_$i.output 2>&1 &
  done

  for job in `jobs -p`
  do
      wait $job
      echo "  job with pid=$job finished"
  done
  
  # STEP 4
  echo ""
  echo "Step 4: Events"
  echo "" >> powheg.input
  echo "#Step 3: Events" >> powheg.input

  overwrite_powheg_var "numevts" $NEVENTS
  overwrite_powheg_var "parallelstage" 4

  echo "  starting $JOBS job(s)..."
  for i in `seq 1 $JOBS`; do
     NSEED=$((i+NSEEDOFFSET))
     echo "  job $i with nseed $NSEED"
     nohup nice -n $NICENESS $EXEPATH <<< $NSEED > $RUNDIR/powheg_step4_$i.output 2>&1 &
  done

  for job in `jobs -p`
  do
      wait $job
      echo "  job with pid=$job finished"
  done

  # experimental
  cat $RUNDIR/pwgevents-*.lhe | grep -v "/LesHouchesEvents" > $RUNDIR/pwgevents.lhe
  echo "</LesHouchesEvents>" >> $RUNDIR/pwgevents.lhe
fi
