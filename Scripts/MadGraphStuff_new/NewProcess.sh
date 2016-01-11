#!/bin/bash

if [[  ! -d ./Cards  ]]; then
    cd ../
    if [[ ! -d ./Cards ]]; then
	echo "Error: cannot find Cards directory"
	exit
    fi
fi

dir=`pwd`
dirm=$dir/MadGraph_POWHEG
main=$dirm/my_proc
rpath=./MadGraph_POWHEG/my_proc

if [[  ! -d $dirm/Template ]]; then
    echo "Error: cannot find MadGraph_POWHEG directory"
    exit
fi
if [[  -d $main ]]; then
    echo 'removing old process directory'
    rm -rf $main
fi

cp -rf $dirm/Template/ $main  >& /dev/null
cp -f $dir/Cards/proc_card.dat $main/Cards >& /dev/null

cd $main
./bin/gen_process $*

cd $dir

cp $rpath/Source/MODEL/coupl.inc . >& /dev/null
cp $rpath/SubProcesses/nlegborn.h . >& /dev/null
cp $rpath/SubProcesses/init_processes.f . >& /dev/null
cp $rpath/lib/libdhelas3.a . >& /dev/null
cp $rpath/lib/libmadgraph.a . >& /dev/null
cp $rpath/lib/libmodel.a . >& /dev/null
cp ../$rpath/Cards/param_card.dat ./Cards/ >& /dev/null
