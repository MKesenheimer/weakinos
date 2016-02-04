#!/bin/bash
#
# Examples:
# runparallel.sh -d testrun_1 -e pwhg_main_nixj
# -> runs pwhg_main_nixj in testrun_1 on 4 cores (default)
#
# runparallel.sh -c -d testrun_1 -e pwhg_main_nixj --itmx1 4 \
# --itmx2 4 --itmx1osres 6 --itmx2osres 8 --ncall1 2000 --ncall2 2000 \
# --ncall1osres 20000 --ncall2osres 20000
# -> runs pwhg_main_nixj in testrun_1 on 4 cores and overwrites some powheg
#    parameters in powheg.input
#
# Hint: use softpoint (must be installed separately) to generate a slha 
#       input file which is then processed with powheg
# softpoint.x sugra --m0=125 --m12=200 --a0=-300 --tanBeta=10 > ./testrun_clean/input.slha 
# ./runparallel.sh -g -c -e pwhg_main_nixj --lopdf 10042 --slha input.slha -d testrun_1
# -> copies the folder testrun_clean (and renames it to testrun_1) and proceeds to
#    calculate the LO cross section on 4 cores.
#
# automatically runs POWHEG in parallel mode with the following stages:
#Stage 1a:
# - fakevirtuals 1
# - xgriditeration 1
# - parallelstage 1
#
#Stage 1b:
# - fakevirtuals 1
# - xgriditeration 2
# - parallelstage 1
#
#Stage 2: 
# - fakevirtuals 0
# - xgriditeration 1
# - parallelstage 2
# - numevts 0
# - nubound 0
#
#Additionally:
#Stage 3:
# - nubound 100000
# - parallelstage 3
#
#Stage 4:
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
  
Optional arguments:
  -h, --help               print this help message
  -i, --info               show informations
  -p, --parallel <n>       number of parallel jobs to submit
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
  -g, --genfolder          generate a new run directory with default input files
  -s, --slha <name>        name of the slha file you want to use
  --lopdf <n>              only LO calculation with LO pdf and LHA number n
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
JOBS=4
CLEAN=false
NSEEDOFFSET=0
NICENESS=10
USEMSUB=false
GENFOLGDER=false

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
    -s|--slha)
        SLHA="$2"
        shift
        shift
        ;;
    --lopdf)
        LOPDF="$2"
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
    --usemsub)
        USEMSUB=true
        shift
        ;;
    -g|--genfolder)
        GENFOLGDER=true
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

# throw a warning if the user wants to use msub
if [ "$USEMSUB" = true ]; then
   echo "Note: use 'msub -I -V -l nodes=1:ppn=$JOBS,walltime=02:00:00' to start an interactive shell and start again."
   echo "Use 'exit' when your calculation is completed."
   echo "Script stopped."
   exit 0
fi

# check if RUNDIR is set
if [ "$RUNDIR" == "" ]; then
   echo "Error: no directory specified."
   usage
fi

#generate a new run directory
if [ "$GENFOLGDER" = true ]; then
   # if the folder already exists
   if [ -d "$RUNDIR" ]; then
      echo "Warning: old directory $RUNDIR will be deleted." 
      read -n1 -r -p "Type [y] to continue... " KEY
      if [ "$KEY" = 'y' ]; then
         rm -r ./$RUNDIR
      else
         echo "Script stopped."
         exit 0
      fi  
   fi
   if [ ! -d "testrun_clean" ]; then
      echo ""
      echo "Error: directory 'testrun_clean' does not exist."
      exit 0
   fi
   cp -r testrun_clean/ ./$RUNDIR
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
if [ "$CLEAN" = true ]; then
   find $RUNDIR ! \( -name '*.slha' -o -name '*.input' -o -name 'pwgseeds.dat' \) -type f -exec rm -f {} +
   cp powheg_clean.input powheg.input
fi
#exit 0

# append to powheg.input
echo "" >> powheg.input
echo "" >> powheg.input
echo "# Modified by runparallel.sh:" >> powheg.input

#default parameters in powheg.input
overwrite_powheg_var "use-old-grid" 1
overwrite_powheg_var "use-old-ubound" 1
comment_powheg_var "iseed"
overwrite_powheg_var "manyseeds" 1
overwrite_powheg_var "testplots" 1

# sed magic
if [ "$NCALL1" != "" ]; then
   overwrite_powheg_var "ncall1" $NCALL1
fi  

if [ "$NCALL2" != "" ]; then
   overwrite_powheg_var "ncall2" $NCALL2
fi

if [ "$NCALL1OSRES" != "" ]; then
   overwrite_powheg_var "ncall1osres" $NCALL1OSRES
fi  

if [ "$NCALL2OSRES" != "" ]; then
   overwrite_powheg_var "ncall2osres" $NCALL2OSRES
fi

if [ "$ITMX1" != "" ]; then
   overwrite_powheg_var "itmx1" $ITMX1
fi  

if [ "$ITMX2" != "" ]; then
   overwrite_powheg_var "itmx2" $ITMX2
fi

if [ "$ITMX1OSRES" != "" ]; then
   overwrite_powheg_var "itmx1osres" $ITMX1OSRES
fi  

if [ "$ITMX2OSRES" != "" ]; then
   overwrite_powheg_var "itmx2osres" $ITMX2OSRES
fi

if [ "$SLHA" != "" ]; then
   overwrite_powheg_var "SLHA" \'"$SLHA"\'
fi

if [ "$LOPDF" != "" ]; then
   overwrite_powheg_var "lhans1" $LOPDF
   overwrite_powheg_var "lhans2" $LOPDF
   overwrite_powheg_var "bornonly" 1
fi

# start the POWHEG-main executable

# STEP 1a
echo ""
echo "Stage 1a: Generating Grids, iteration 1"
echo "" >> powheg.input
echo "#Stage 1a: Generating Grids, iteration 1" >> powheg.input

overwrite_powheg_var "fakevirtuals" 1
overwrite_powheg_var "parallelstage" 1
overwrite_powheg_var "xgriditeration" 1

echo "  starting $JOBS job(s)..."
for i in `seq 1 $JOBS`; do
   NSEED=$((i+NSEEDOFFSET))
   echo "  job $i with nseed $NSEED"
   nohup nice -n $NICENESS $EXEPATH <<< $NSEED > $RUNDIR/powheg_step1a_$i.output 2>&1 &
done

for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
done

# STEP 1b
echo ""
echo "Stage 1b: Generating Grids, iteration 2"
echo "" >> powheg.input
echo "#Stage 1b: Generating Grids, iteration 2" >> powheg.input

overwrite_powheg_var "fakevirtuals" 1
overwrite_powheg_var "parallelstage" 1
overwrite_powheg_var "xgriditeration" 2

echo "  starting $JOBS job(s)..."
for i in `seq 1 $JOBS`; do
   NSEED=$((i+NSEEDOFFSET))
   echo "  job $i with nseed $NSEED"
   nohup nice -n $NICENESS $EXEPATH <<< $NSEED > $RUNDIR/powheg_step1b_$i.output 2>&1 &
done

for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
done

# STEP 2
echo ""
echo "Stage 2: NLO run"
echo "" >> powheg.input
echo "#Stage 2: NLO run" >> powheg.input

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

for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
done

# combined results for stage 2
echo ""
echo "Combined results for stage 2:"
#$(find \( -name 'pwg-*-stat.dat' -and ! \( -name 'pwg-st2*' \) \) -type f | awk '{print}' ORS=' ')
$WORKINGDIR/merge-pwg-stat $(ls ./pwg-st2-*-stat.dat) > $RUNDIR/pwg-st2-combined-stat.dat
cat $RUNDIR/pwg-st2-combined-stat.dat

# if the user wants to generate events
if [ "$NEVENTS" != "" ]; then
  # STEP 3
  echo ""
  echo "Stage 3: Upper bound"
  echo "" >> powheg.input
  echo "#Stage 3: Upper bound" >> powheg.input

  overwrite_powheg_var "nubound" $NUBOUND
  overwrite_powheg_var "parallelstage" 3

  echo "  starting $JOBS job(s)..."
  for i in `seq 1 $JOBS`; do
     NSEED=$((i+NSEEDOFFSET))
     echo "  job $i with nseed $NSEED"
     nohup nice -n $NICENESS $EXEPATH <<< $NSEED > $RUNDIR/powheg_step3_$i.output 2>&1 &
  done

  for job in `jobs -p`; do
      wait $job
      echo "  job with pid=$job finished"
  done
  
  # STEP 4
  echo ""
  echo "Stage 4: Events"
  echo "" >> powheg.input
  echo "#Stage 4: Events" >> powheg.input

  overwrite_powheg_var "numevts" $NEVENTS
  overwrite_powheg_var "parallelstage" 4

  echo "  starting $JOBS job(s)..."
  for i in `seq 1 $JOBS`; do
     NSEED=$((i+NSEEDOFFSET))
     echo "  job $i with nseed $NSEED"
     nohup nice -n $NICENESS $EXEPATH <<< $NSEED > $RUNDIR/powheg_step4_$i.output 2>&1 &
  done

  for job in `jobs -p`; do
      wait $job
      echo "  job with pid=$job finished"
  done

  # experimental
  cat $RUNDIR/pwgevents-*.lhe | grep -v "/LesHouchesEvents" > $RUNDIR/pwgevents.lhe
  echo "</LesHouchesEvents>" >> $RUNDIR/pwgevents.lhe
fi
