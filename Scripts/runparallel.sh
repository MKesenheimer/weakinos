#!/bin/bash
#
# Examples:
#
# $ ./runparallel.sh -d testrun_1 -e pwhg_main_nixj
# -> runs pwhg_main_nixj in testrun_1 on 4 cores (default)
#
#
# $ ./runparallel.sh -c -d testrun_1 -e pwhg_main_nixj --itmx1 4 \
# --itmx2 4 --itmx1osres 6 --itmx2osres 8 --ncall1 2000 --ncall2 2000 \
# --ncall1osres 20000 --ncall2osres 20000
# -> runs pwhg_main_nixj in testrun_1 on 4 cores and overwrites some powheg
#    parameters in powheg.input
#
#
# $ ./runparallel.sh -g -c -e pwhg_main_nixj -d run_wevents --genevents > log_wevents
# -> copies the folder testrun_clean (and renames it to run_wevents),
#    generates events (nubound and nevents in powheg.input must be greater than zero,
#    or use --nevents and --nubound to set the numbers)
#
#
# $ softpoint.x sugra --m0=125 --m12=200 --a0=-300 --tanBeta=10 > ./testrun_clean/input.slha 
# $ ./runparallel.sh -g -c -e pwhg_main_nixj --lopdf 10042 --slha input.slha -d testrun_1
# -> use softpoint (must be installed separately) to generate a slha input file which
#    is then processed with powheg.
#    Copies the folder testrun_clean (and renames it to testrun_1) and proceeds to
#    calculate the LO cross section on 4 cores.


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
  --genevents              generate events
  --nubound <n>            upper bound number
  --nevents <n>            number of events
  --usemsub                use the submitting system msub (not implemented yet)
  -s, --slha <name>        name of the slha file you want to use
  --lopdf <n>              only LO calculation with LO pdf and LHA number n
  -g, --genfolder          generate a new run directory with default input files 
                           (the directory "testrun_clean" is needed)
  --fin1 <n>               PDG number of first final particle
  --fin2 <n>               PDG number of second final particle
  --dec1 <n>               PDG number of first decay particle
  --dec2 <n>               PDG number of second decay particle
  --dec3 <n>               PDG number of third decay particle
  --dec4 <n>               PDG number of fourth decay particle
  --mur <n>                set renscfact
  --muf <n>                set facscfact
  --ewi <n>                set the regulator ewi
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

function checksum {
  if type md5 1>/dev/null 2>/dev/null; then
    md5 "$@"
  elif type md5sum 1>/dev/null 2>/dev/null; then
    md5sum "$@"
  elif type sha1sum 1>/dev/null 2>/dev/null; then
    sha1sum "$@"
  elif type sha256sum 1>/dev/null 2>/dev/null; then
    sha256sum "$@"
  else
    echo $(date +%s)
  fi
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
GENEVENTS=false
GENFOLGDER=false
MERGE=false
ARG1=""

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
        GENEVENTS=true
        shift
        ;;
    --nubound)
        GENEVENTS=true
        NUBOUND="$2"
        shift
        shift
        ;;
    --nevents)
        GENEVENTS=true
        NEVENTS="$2"
        shift
        shift
        ;;
    --fin1)
        FIN1="$2"
        shift
        shift
        ;;
    --fin2)
        FIN2="$2"
        shift
        shift
        ;;
    --dec1)
        DEC1="$2"
        shift
        shift
        ;;
    --dec2)
        DEC2="$2"
        shift
        shift
        ;;
    --dec3)
        DEC3="$2"
        shift
        shift
        ;;
    --dec4)
        DEC4="$2"
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
    --ewi)
        EWI="$2"
        shift
        shift
        ;;
    --usemsub)
        USEMSUB=true
        shift
        ;;
    --merge)
        MERGE=true
        shift
        ;;
    -g|--genfolder)
        GENFOLGDER=true
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

# directories
WORKINGDIR=${PWD}
RUNDIR=$WORKINGDIR/$RUNDIR
EXEPATH=$WORKINGDIR/$EXE

#generate a new run directory
if [ "$GENFOLGDER" = true ]; then
   # if the folder already exists
   if [ -d "$RUNDIR" ]; then
      echo "Warning: old directory $RUNDIR will be deleted." 
      read -n1 -r -p "Type [y] to continue... " KEY
      if [ "$KEY" = 'y' ]; then
         rm -r $RUNDIR
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
   cp -r testrun_clean/ $RUNDIR
   cp $RUNDIR/powheg_clean.input $RUNDIR/powheg.input
fi

# check if EXE is set
if [ "$EXE" == "" ]
then
   echo "Error: no executable specified."
   usage
fi

# generate an individual identifier for every run
IDENT=$(echo -n "$RUNDIR" | checksum | cut -c1-8)
echo ""
echo "Identifier is $IDENT"

# change into the directory where to execute pwhg_main
cd $RUNDIR

# clean up the directory
if [ "$CLEAN" = true ]; then
   find $RUNDIR ! \( -name '*.slha' -o -name '*.input' -o -name 'pwgseeds.dat' \) -type f -exec rm -f {} +
   cp $RUNDIR/powheg_clean.input $RUNDIR/powheg.input
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
   overwrite_powheg_var "LOevents" 1
fi

if [ "$FIN1" != "" ]; then
   overwrite_powheg_var "fin1" "$FIN1"
fi

if [ "$FIN2" != "" ]; then
   overwrite_powheg_var "fin2" "$FIN2"
fi

if [ "$DEC1" != "" ]; then
   overwrite_powheg_var "dec1" "$DEC1"
fi

if [ "$DEC2" != "" ]; then
   overwrite_powheg_var "dec2" "$DEC2"
fi

if [ "$DEC3" != "" ]; then
   overwrite_powheg_var "dec3" "$DEC3"
fi

if [ "$DEC4" != "" ]; then
   overwrite_powheg_var "dec4" "$DEC4"
fi

if [ "$RENSCFACT" != "" ]; then
   overwrite_powheg_var "renscfact" "$RENSCFACT"
fi

if [ "$FACSCFACT" != "" ]; then
   overwrite_powheg_var "facscfact" "$FACSCFACT"
fi

if [ "$EWI" != "" ]; then
   overwrite_powheg_var "ewi" "$EWI"
fi

# back up the old parameters
NEVENTSOLD=$(read_powheg_var "numevts")
NUBOUNDOLD=$(read_powheg_var "nubound")

#echo "Note: numevts was $NEVENTSOLD"
#echo "Note: nubound was $NUBOUNDOLD"

# generate the scripts to start the POWHEG-main executable

# STEP 1a
echo "" >> $RUNDIR/powheg.input
echo "#Stage 1a: Generating Grids, iteration 1" >> $RUNDIR/powheg.input
overwrite_powheg_var "fakevirtuals" 1
overwrite_powheg_var "parallelstage" 1
overwrite_powheg_var "xgriditeration" 1
cp $RUNDIR/powheg.input $RUNDIR/powheg_st1a.input

# wheter \$1 or \$ARG1 is defined (msub sets ARG1)
cat <<EOM > $WORKINGDIR/run_st1a_${IDENT}.sh
#!/bin/bash
cd $RUNDIR
cp $RUNDIR/powheg_st1a.input $RUNDIR/powheg.input
$EXEPATH < <(printf "%s\n" "\$1" "\$ARG1")
EOM
chmod +x $WORKINGDIR/run_st1a_${IDENT}.sh

# STEP 1b
echo "" >> $RUNDIR/powheg.input
echo "#Stage 1b: Generating Grids, iteration 2" >> $RUNDIR/powheg.input
overwrite_powheg_var "fakevirtuals" 1
overwrite_powheg_var "parallelstage" 1
overwrite_powheg_var "xgriditeration" 2
cp $RUNDIR/powheg.input $RUNDIR/powheg_st1b.input

cat <<EOM > $WORKINGDIR/run_st1b_${IDENT}.sh
#!/bin/bash
cd $RUNDIR
cp $RUNDIR/powheg_st1b.input $RUNDIR/powheg.input
$EXEPATH < <(printf "%s\n" "\$1" "\$ARG1")
EOM
chmod +x $WORKINGDIR/run_st1b_${IDENT}.sh

# STEP 2
echo "" >> $RUNDIR/powheg.input
echo "#Stage 2: NLO run" >> $RUNDIR/powheg.input
overwrite_powheg_var "fakevirtuals" 0
# TODO: xgriditeration = 1 necessary?
overwrite_powheg_var "xgriditeration" 1
overwrite_powheg_var "parallelstage" 2
overwrite_powheg_var "numevts" 0
overwrite_powheg_var "nubound" 0
cp $RUNDIR/powheg.input $RUNDIR/powheg_st2.input

cat <<EOM > $WORKINGDIR/run_st2_${IDENT}.sh
#!/bin/bash
cd $RUNDIR
cp $RUNDIR/powheg_st2.input $RUNDIR/powheg.input
$EXEPATH < <(printf "%s\n" "\$1" "\$ARG1")
EOM
chmod +x $WORKINGDIR/run_st2_${IDENT}.sh

# if the user wants to generate events
if [ "$GENEVENTS" = true ]; then

# STEP 3
echo "" >> $RUNDIR/powheg.input
echo "#Stage 3: Upper bound" >> $RUNDIR/powheg.input
if [ "$NUBOUND" != "" ]; then
  overwrite_powheg_var "nubound" $NUBOUND
else
  overwrite_powheg_var "nubound" $NUBOUNDOLD
fi
overwrite_powheg_var "parallelstage" 3
cp $RUNDIR/powheg.input $RUNDIR/powheg_st3.input

cat <<EOM > $WORKINGDIR/run_st3_${IDENT}.sh
#!/bin/bash
cd $RUNDIR
cp $RUNDIR/powheg_st3.input $RUNDIR/powheg.input
$EXEPATH < <(printf "%s\n" "\$1" "\$ARG1")
EOM
chmod +x $WORKINGDIR/run_st3_${IDENT}.sh

# STEP 4
echo "" >> $RUNDIR/powheg.input
echo "#Stage 4: Events" >> $RUNDIR/powheg.input
if [ "$NEVENTS" != "" ]; then
  overwrite_powheg_var "numevts" $NEVENTS
else
  overwrite_powheg_var "numevts" $NEVENTSOLD
fi
overwrite_powheg_var "parallelstage" 4
cp $RUNDIR/powheg.input $RUNDIR/powheg_st4.input

cat <<EOM > $WORKINGDIR/run_st4_${IDENT}.sh
#!/bin/bash
cd $RUNDIR
cp $RUNDIR/powheg_st4.input $RUNDIR/powheg.input
$EXEPATH < <(printf "%s\n" "\$1" "\$ARG1")
EOM
chmod +x $WORKINGDIR/run_st4_${IDENT}.sh
fi

# generate and run the run.sh script
# nohup ./runparallel.sh -g -c -e pwhg_main_nixj -d run_nsusy_n2x1+ -p 4 --fin1 1000023 --fin2 1000024 --slha input_nsusy_1307.0782.slha --ncall1 20000 --ncall2 20000 --ncall1osres 2000000 --ncall2osres 2000000 --nevents 500000 --nubound 500000 --genevents --merge > log_run1_nsusy_n2x1+ &
# nohup ./runparallel.sh -g -c -e pwhg_main_nixj -d run_mSUGRA_n2x1+ -p 4 --fin1 1000023 --fin2 1000024 --slha input_mSUGRA_1410.4999.slha --ncall1 20000 --ncall2 20000 --ncall1osres 2000000 --ncall2osres 2000000 --nevents 500000 --nubound 500000 --genevents --merge > log_run_mSUGRA_n2x1+ &
if [ "$USEMSUB" = false ]; then
cat <<EOM > $WORKINGDIR/run_${IDENT}.sh
#!/bin/bash
echo ""
echo "Stage 1a: Generating Grids, iteration 1"
echo "  starting $JOBS job(s)..."
for i in \`seq 1 $JOBS\`; do
  NSEED=\$((\$i+$NSEEDOFFSET))
  echo "  job \$i with nseed \$NSEED"
  nohup nice -n $NICENESS $WORKINGDIR/run_st1a_${IDENT}.sh \$NSEED > $RUNDIR/powheg_st1a_\${NSEED}.output 2>&1 &
done
for job in \`jobs -p\`; do
    wait \$job
    echo "  job with pid=\$job finished"
done

echo ""
echo "Stage 1b: Generating Grids, iteration 2"
echo "  starting $JOBS job(s)..."
for i in \`seq 1 $JOBS\`; do
   NSEED=\$((i+$NSEEDOFFSET))
   echo "  job \$i with nseed \$NSEED"
   nohup nice -n $NICENESS $WORKINGDIR/run_st1b_${IDENT}.sh \$NSEED > $RUNDIR/powheg_st1b_\${NSEED}.output 2>&1 &
done
for job in \`jobs -p\`; do
    wait \$job
    echo "  job with pid=\$job finished"
done

echo ""
echo "Stage 2: NLO run"
echo "  starting $JOBS job(s)..."
for i in \`seq 1 $JOBS\`; do
   NSEED=\$((i+$NSEEDOFFSET))
   echo "  job \$i with nseed \$NSEED"
   nohup nice -n $NICENESS $WORKINGDIR/run_st2_${IDENT}.sh \$NSEED > $RUNDIR/powheg_st2_\${NSEED}.output 2>&1 &
done
for job in \`jobs -p\`; do
    wait \$job
    echo "  job with pid=\$job finished"
done
# combined results for stage 2
echo ""
echo "Combined results for stage 2:"
rm $RUNDIR/pwg-st2-combined-stat.dat
cd $RUNDIR && ../merge-pwg-stat \$(ls ./pwg-st2-*-stat.dat) > pwg-st2-combined-stat.dat
cat $RUNDIR/pwg-st2-combined-stat.dat
EOM

if [ "$GENEVENTS" = true ]; then
cat <<EOM >> $WORKINGDIR/run_${IDENT}.sh
echo ""
echo "Stage 3: Upper bound"
echo "  starting $JOBS job(s)..."
for i in \`seq 1 $JOBS\`; do
   NSEED=\$((i+$NSEEDOFFSET))
   echo "  job \$i with nseed \$NSEED"
   nohup nice -n $NICENESS $WORKINGDIR/run_st3_${IDENT}.sh \$NSEED > $RUNDIR/powheg_st3_\${NSEED}.output 2>&1 &
done
for job in \`jobs -p\`; do
    wait \$job
    echo "  job with pid=\$job finished"
done

echo ""
echo "Stage 4: Events"
echo "  starting $JOBS job(s)..."
for i in \`seq 1 $JOBS\`; do
   NSEED=\$((i+$NSEEDOFFSET))
   echo "  job \$i with nseed \$NSEED"
   nohup nice -n $NICENESS $WORKINGDIR/run_st4_${IDENT}.sh \$NSEED > $RUNDIR/powheg_st4_\${NSEED}.output 2>&1 &
done
for job in \`jobs -p\`; do
    wait \$job
    echo "  job with pid=\$job finished"
done

EOM

if [ "$MERGE" = true ]; then
cat <<EOM >> $WORKINGDIR/run_${IDENT}.sh
# merge the event files
cat $RUNDIR/pwgevents-*.lhe | grep -v "/LesHouchesEvents" > $RUNDIR/pwgevents.lhe
echo "</LesHouchesEvents>" >> $RUNDIR/pwgevents.lhe
#if [ -e "$RUNDIR/pwgevents.lhe" ]; then
#  echo "merged event files succesfully, deleting old event files..."
#  find $RUNDIR -type f -name "pwgevents-*" -exec rm -f '{}' \;
#fi
# merge the NLO top files
rm $RUNDIR/pwg-NLO.top
cd $RUNDIR && ../merge-data 1 \$(ls ./pwg-*-NLO.top) && mv fort.12 pwg-NLO.top

EOM
fi #if MERGE
fi #if GENEVENTS

chmod +x $WORKINGDIR/run_${IDENT}.sh
$WORKINGDIR/run_${IDENT}.sh

# if finished delete the old files
rm $WORKINGDIR/run_st1a_${IDENT}.sh
rm $WORKINGDIR/run_st1b_${IDENT}.sh
rm $WORKINGDIR/run_st2_${IDENT}.sh
rm $WORKINGDIR/run_st3_${IDENT}.sh
rm $WORKINGDIR/run_st4_${IDENT}.sh
rm $WORKINGDIR/run_${IDENT}.sh
rm $RUNDIR/powheg_st*.input
fi


# if the user wants to use msub:
# the approximate runtime is determined with the following parameters for the nemo cluster in freiburg:
# low precision job: 20 parallel jobs
# ./runparallel.sh -g -c -e pwhg_main_nixj -d run_nsusy_n2x1+ -p 20 --fin1 1000023 --fin2 1000024 --slha input_nsusy_1307.0782.slha --nevents 50000 --nubound 50000 --genevents --usemsub
# ncall1   20000
# itmx1    4
# ncall2   20000
# itmx2    4
# ncall1osres 2000000
# itmx1osres  6
# ncall2osres 2000000
# itmx2osres  8 
# nubound 50000
# numevts 50000
#
# stage 1a: 8-12min
# stage 1b: 8-12min
# stage 2: 2-3h
# stage 3: 20-30sec
# stage 4: 15-25min
# total: ~4h

# high precision job: 20 parallel jobs
# ./runparallel.sh -g -c -e pwhg_main_nixj -d run1_nsusy_n2x1+ -p 20 --fin1 1000023 --fin2 1000024 --slha input_nsusy_1307.0782.slha --ncall1 200000 --ncall2 300000 --nevents 100000 --nubound 100000 --genevents --usemsub > submit_run1_nsusy_n2x1+ &
# ./runparallel.sh -g -c -e pwhg_main_nixj -d run2_nsusy_n2x1+ -p 20 --fin1 1000023 --fin2 1000024 --slha input_nsusy_1307.0782.slha --ncall1 200000 --ncall2 300000 --nevents 100000 --nubound 100000 --genevents --usemsub --offset 20 > submit_run2_nsusy_n2x1+ &
# tail -f submit_run_old_nsusy_n2x1+
# ncall1 200000
# itmx1    4
# ncall2 300000
# itmx2    4
# ncall1osres 2000000
# itmx1osres  6
# ncall2osres 2000000
# itmx2osres  8 
# nubound 100000
# numevts 100000
#
# stage 1a: 10min
# stage 1b: 15min
# stage 2: 5h
# stage 3: 3min
# stage 4: 40min
# total: ~6h

if [ "$USEMSUB" = true ]; then
cat <<EOM > $WORKINGDIR/runmsub_${IDENT}.sh
#!/bin/bash
echo ""
echo "Stage 1a: Generating Grids, iteration 1"
echo "  submitting $JOBS job(s)..."
for i in \`seq 1 $JOBS\`; do
  NSEED=\$((\$i+$NSEEDOFFSET))
  job[\$i]=\$(msub -l walltime=01:00:00 -v ARG1=\$NSEED -o $RUNDIR/powheg_st1a_\${NSEED}.output -e $RUNDIR/powheg_st1a_\${NSEED}.error $WORKINGDIR/run_st1a_${IDENT}.sh | grep -v -e '^$')
  echo "  job \$i with nseed \$NSEED and ID \${job[\$i]}"
  dependIDs1a="\$dependIDs1a:\${job[\$i]}"
  #echo \$dependIDs
done

echo ""
echo "Stage 1b: Generating Grids, iteration 2"
echo "  submitting $JOBS job(s)..."
for i in \`seq 1 $JOBS\`; do
  NSEED=\$((\$i+$NSEEDOFFSET))
  job[\$i]=\$(msub -l walltime=01:00:00,depend=afterok\${dependIDs1a} -v ARG1=\$NSEED -o $RUNDIR/powheg_st1b_\${NSEED}.output -e $RUNDIR/powheg_st1b_\${NSEED}.error $WORKINGDIR/run_st1b_${IDENT}.sh | grep -v -e '^$')
  echo "  job \$i with nseed \$NSEED and ID \${job[\$i]}"
  dependIDs1b="\$dependIDs1b:\${job[\$i]}"
done

echo ""
echo "Stage 2: NLO run"
echo "  submitting $JOBS job(s)..."
for i in \`seq 1 $JOBS\`; do
  NSEED=\$((\$i+$NSEEDOFFSET))
  job[\$i]=\$(msub -l walltime=12:00:00,depend=afterok\${dependIDs1b} -v ARG1=\$NSEED -o $RUNDIR/powheg_st2_\${NSEED}.output -e $RUNDIR/powheg_st2_\${NSEED}.error $WORKINGDIR/run_st2_${IDENT}.sh | grep -v -e '^$')
  echo "  job \$i with nseed \$NSEED and ID \${job[\$i]}"
  dependIDs2="\$dependIDs2:\${job[\$i]}"
done

EOM
if [ "$GENEVENTS" = true ]; then
cat <<EOM >> $WORKINGDIR/runmsub_${IDENT}.sh
echo ""
echo "Stage 3: Upper bound"
echo "  submitting $JOBS job(s)..."
for i in \`seq 1 $JOBS\`; do
  NSEED=\$((\$i+$NSEEDOFFSET))
  job[\$i]=\$(msub -l walltime=00:30:00,depend=afterok\${dependIDs2} -v ARG1=\$NSEED -o $RUNDIR/powheg_st3_\${NSEED}.output -e $RUNDIR/powheg_st3_\${NSEED}.error $WORKINGDIR/run_st3_${IDENT}.sh | grep -v -e '^$')
  echo "  job \$i with nseed \$NSEED and ID \${job[\$i]}"
  dependIDs3="\$dependIDs3:\${job[\$i]}"
done

echo ""
echo "Stage 4: Events"
echo "  submitting $JOBS job(s)..."
for i in \`seq 1 $JOBS\`; do
  NSEED=\$((\$i+$NSEEDOFFSET))
  job[\$i]=\$(msub -l walltime=12:00:00,depend=afterok\${dependIDs3} -v ARG1=\$NSEED -o $RUNDIR/powheg_st4_\${NSEED}.output -e $RUNDIR/powheg_st4_\${NSEED}.error $WORKINGDIR/run_st4_${IDENT}.sh | grep -v -e '^$')
  echo "  job \$i with nseed \$NSEED and ID \${job[\$i]}"
  dependIDs4="\$dependIDs4:\${job[\$i]}"
done

EOM
fi #if GENEVENTS

chmod +x $WORKINGDIR/runmsub_${IDENT}.sh
$WORKINGDIR/runmsub_${IDENT}.sh
fi