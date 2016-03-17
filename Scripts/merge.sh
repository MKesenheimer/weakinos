#!/bin/bash
# Usage: ./merge.sh <directory>
WORKINGDIR=${PWD}
RUNDIR=$WORKINGDIR/$1

# combined results for stage 2
echo ""
echo "Combined results for stage 2:"
$WORKINGDIR/merge-pwg-stat $(ls $RUNDIR/pwg-st2-*-stat.dat) > $RUNDIR/pwg-st2-combined-stat.dat
cat $RUNDIR/pwg-st2-combined-stat.dat

# merge the event files
echo ""
echo "Merging event files..."
cat $RUNDIR/pwgevents-*.lhe | grep -v "/LesHouchesEvents" > $RUNDIR/pwgevents.lhe
echo "</LesHouchesEvents>" >> $RUNDIR/pwgevents.lhe
#if [ -e "$RUNDIR/pwgevents.lhe" ]; then
#  echo "merged event files succesfully, deleting old event files..."
#  find $RUNDIR -type f -name "pwgevents-*" -exec rm -f '{}' \;
#fi
# merge the NLO top files
cd $RUNDIR && ../merge-data 1 $(ls $RUNDIR/pwg-*-NLO.top) && mv fort.12 pwg-NLO.top
echo "done."