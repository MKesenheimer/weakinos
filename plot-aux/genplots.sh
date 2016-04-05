#!/bin/bash

function exitmessage {
echo " genplot called with $narg arguments $fil1 $fil2 $out"
echo "usage: genplots.sh file1 file2 nameoutput"
echo "   or: genplots.sh file nameoutput"
exit 0
}

function check {
if ! [ -e $1 ]
then
    echo file $1 not found
    exitmessage
fi
}

case $# in
 2) fil2=$1 ; out=$2 ; check $fil2 ;;
 3) fil1=$1 ; fil2=$2 ; out=$3 ; check $fil1 ; check $fil2 ;;
 *) exitmessage ;;
esac

>genplots.gp
i=0
cat $fil2 | grep '#' | while read line
do
title=`echo $line |sed 's/#// ; s/ *index *[0-9]* *// ; s/ *//g'`
echo "$title"

ttt=`echo "$line"|sed 's/#//g ; s/ /_/g ; s/\//o/g ' `

file=$out-$ttt

cat <<EOF >> genplots.gp
reset
set term post eps enhanced color dashed "Helvetica" 24
set output '$file.eps'
iindex = $i
mytitle = "$title"
EOF
if ! [ a$fil1 = a ]
then
cat <<EOF >> genplots.gp
infil1 = "$fil1"
infil2 = "$fil2"
#both = '< paste $fil1 $fil2'
both = "< ../../plot-aux/pastegnudata "."'"."[$title]"."'"." $fil1 $fil2"
set multiplot

set origin 0.0,0.3
set size 1.0,0.7
set lmargin at screen 0.2
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.4
EOF
else
cat <<EOF >> genplots.gp
infil1 = "$fil2"
EOF
fi
cat <<EOF >> genplots.gp

#set log x
set log y

set datafile fortran
set style data xyerrorlines

set xtics format ""

set xrange [] writeback
set yrange [] 
set ylabel '{d {/Symbol s}}/{d {p_T}} [fb]'
set format y "10^{%L}"

EOF

if [ a$fil1 = a ]
then
cat <<'EOF' >> genplots.gp
plot infil1 index  iindex using ($1+$2)/2:3:1:2:($3-$4):($3+$4) title mytitle

EOF
else
cat <<'EOF' >> genplots.gp
plot infil1 index  iindex using ($1+$2)/2:(1*$3):1:2:(1*($3-$4)):(1*($3+$4)) title mytitle,\
    infil2 index  iindex  using ($1+$2)/2:(1*$3):1:2:(1*($3-$4)):(1*($3+$4)) title mytitle

set tmargin at screen 0.4
set bmargin at screen 0.15

set nolog y
set yrange [0:2]
set xrange restore
set format x
set format y
set xlabel '[GeV]'
set ylabel ''
set ytics 0.5,0.5,1.5
plot both index 0 \
using ($1+$2)/2:($3/$7):1:2:($3/$7)*(1-(($4/$3)**2+($8/$7)**2)**0.5):($3/$7)*(1+(($4/$3)**2+($8/$7)**2)**0.5) title '', 1 with lines linestyle -1 title ''
unset key

unset multiplot

EOF
fi
i=$[$i+1]
done
