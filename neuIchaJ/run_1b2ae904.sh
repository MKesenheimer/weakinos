#!/bin/bash
echo ""
echo "Stage 1a: Generating Grids, iteration 1"
echo "  starting 4 job(s)..."
for i in `seq 1 4`; do
  NSEED=$(($i+0))
  echo "  job $i with nseed $NSEED"
  nohup nice -n 10 /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/run_st1a_1b2ae904.sh $NSEED > /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1/powheg_st1a_${NSEED}.output 2>&1 &
done
for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
done

echo ""
echo "Stage 1b: Generating Grids, iteration 2"
echo "  starting 4 job(s)..."
for i in `seq 1 4`; do
   NSEED=$((i+0))
   echo "  job $i with nseed $NSEED"
   nohup nice -n 10 /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/run_st1b_1b2ae904.sh $NSEED > /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1/powheg_st1b_${NSEED}.output 2>&1 &
done
for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
done

echo ""
echo "Stage 2: NLO run"
echo "  starting 4 job(s)..."
for i in `seq 1 4`; do
   NSEED=$((i+0))
   echo "  job $i with nseed $NSEED"
   nohup nice -n 10 /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/run_st2_1b2ae904.sh $NSEED > /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1/powheg_st2_${NSEED}.output 2>&1 &
done
for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
done
# combined results for stage 2
echo ""
echo "Combined results for stage 2:"
rm /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1/pwg-st2-combined-stat.dat
cd /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1 && ../merge-pwg-stat $(ls ./pwg-st2-*-stat.dat) > pwg-st2-combined-stat.dat
cat /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1/pwg-st2-combined-stat.dat
echo ""
echo "Stage 3: Upper bound"
echo "  starting 4 job(s)..."
for i in `seq 1 4`; do
   NSEED=$((i+0))
   echo "  job $i with nseed $NSEED"
   nohup nice -n 10 /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/run_st3_1b2ae904.sh $NSEED > /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1/powheg_st3_${NSEED}.output 2>&1 &
done
for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
done

echo ""
echo "Stage 4: Events"
echo "  starting 4 job(s)..."
for i in `seq 1 4`; do
   NSEED=$((i+0))
   echo "  job $i with nseed $NSEED"
   nohup nice -n 10 /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/run_st4_1b2ae904.sh $NSEED > /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1/powheg_st4_${NSEED}.output 2>&1 &
done
for job in `jobs -p`; do
    wait $job
    echo "  job with pid=$job finished"
done

# merge the event files
cat /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1/pwgevents-*.lhe | grep -v "/LesHouchesEvents" > /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1/pwgevents.lhe
echo "</LesHouchesEvents>" >> /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1/pwgevents.lhe
#if [ -e "/home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1/pwgevents.lhe" ]; then
#  echo "merged event files succesfully, deleting old event files..."
#  find /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1 -type f -name "pwgevents-*" -exec rm -f '{}' \;
#fi
# merge the NLO top files
rm /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1/pwg-NLO.top
cd /home/hepuser/kesenheimer/POWHEG-BOX-V2/weakinos/neuIchaJ/testrun1 && ../merge-data 1 $(ls ./pwg-*-NLO.top) && mv fort.12 pwg-NLO.top

