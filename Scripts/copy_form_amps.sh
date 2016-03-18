#!/bin/bash

WORKINGDIR=${PWD}
PROCDIR="neuIchaJ"
INITIAL1="dubar"
INITIAL2="udbar"
FINAL="n1x1_test"

# copies the squaredME files generated with FormCalc to the right places
# and calls the rename scripts to rename the files and their contents
echo "Generating directories..."
mkdir ${WORKINGDIR}/../${PROCDIR}/FormCalc_Virtuals/${INITIAL1}_${FINAL}_squaredme/
mkdir ${WORKINGDIR}/../${PROCDIR}/FormCalc_Virtuals/${INITIAL2}_${FINAL}_squaredme/

echo ""
echo "Deleting old files..."
rm ${WORKINGDIR}/../${PROCDIR}/FormCalc_Virtuals/${INITIAL1}_${FINAL}_squaredme/*.F
rm ${WORKINGDIR}/../${PROCDIR}/FormCalc_Virtuals/${INITIAL2}_${FINAL}_squaredme/*.F

rm ${WORKINGDIR}/../${PROCDIR}/FormCalc_Virtuals/include/${INITIAL1}_${FINAL}_vars.h
rm ${WORKINGDIR}/../${PROCDIR}/FormCalc_Virtuals/include/${INITIAL2}_${FINAL}_vars.h

echo ""
echo "Copying new files..."
cp ${WORKINGDIR}/${INITIAL1}.fortran/squaredme/*.F ../${PROCDIR}/FormCalc_Virtuals/${INITIAL1}_${FINAL}_squaredme/
cp ${WORKINGDIR}/${INITIAL2}.fortran/squaredme/*.F ../${PROCDIR}/FormCalc_Virtuals/${INITIAL2}_${FINAL}_squaredme/

cp ${WORKINGDIR}/${INITIAL1}.fortran/squaredme/vars.h ../${PROCDIR}/FormCalc_Virtuals/include/${INITIAL1}_${FINAL}_vars.h
cp ${WORKINGDIR}/${INITIAL2}.fortran/squaredme/vars.h ../${PROCDIR}/FormCalc_Virtuals/include/${INITIAL2}_${FINAL}_vars.h

# call rename scripts
echo ""
echo "Renaming files..."

echo "cd ${WORKINGDIR}/../${PROCDIR}/FormCalc_Virtuals/${INITIAL1}_${FINAL}_squaredme"
cd ${WORKINGDIR}/../${PROCDIR}/FormCalc_Virtuals/${INITIAL1}_${FINAL}_squaredme
PRE="${INITIAL1}_${FINAL}_"
for file in *.F *.f; do
    # Dateien umbenennen -> Prefix anhängen
    echo "${PRE}${file}"
    mv "$file" "${PRE}${file}"
done
for file in $(find . -type f -iname "*.f" -o -iname "*.F"); do
    # alle Vorkommnisse von vars.h durch xxbar_vars.h ersetzen usw.
    sed -i -e "s/vars.h/${PRE}vars.h/g" $file
    sed -i -e "s/born/${PRE}born/g" $file
    sed -i -e "s/vert/${PRE}vert/g" $file
    sed -i -e "s/self/${PRE}self/g" $file
    sed -i -e "s/box/${PRE}box/g" $file
    sed -i -e "s/SquaredME/${PRE}SquaredME/g" $file
    sed -i -e "s/abbr0h/${PRE}abbr0h/g" $file
    sed -i -e "s/abbr1h/${PRE}abbr1h/g" $file
    sed -i -e "s/abbr0s/${PRE}abbr0s/g" $file
    sed -i -e "s/abbr1s/${PRE}abbr1s/g" $file
    sed -i -e "s/abbr1a/${PRE}abbr1a/g" $file    
    sed -i -e 's/#ifdef DEBUG/#ifdef DEBUGQ/g' $file
    # jede Zeile löschen, die #include "inline.h" enthält
    sed -i -e '/^#include "inline.h"/d' $file
done

echo ""
echo "cd ${WORKINGDIR}/../${PROCDIR}/FormCalc_Virtuals/${INITIAL1}_${FINAL}_squaredme"
cd ${WORKINGDIR}/../${PROCDIR}/FormCalc_Virtuals/${INITIAL2}_${FINAL}_squaredme
PRE="${INITIAL2}_${FINAL}_"
for file in *.F *.f; do
    # Dateien umbenennen -> Prefix anhängen
    echo "${PRE}${file}"
    mv "$file" "${PRE}${file}"
done
for file in $(find . -type f -iname "*.f" -o -iname "*.F"); do
    # alle Vorkommnisse von vars.h durch xxbar_vars.h ersetzen usw.
    sed -i -e "s/vars.h/${PRE}vars.h/g" $file
    sed -i -e "s/born/${PRE}born/g" $file
    sed -i -e "s/vert/${PRE}vert/g" $file
    sed -i -e "s/self/${PRE}self/g" $file
    sed -i -e "s/box/${PRE}box/g" $file
    sed -i -e "s/SquaredME/${PRE}SquaredME/g" $file
    sed -i -e "s/abbr0h/${PRE}abbr0h/g" $file
    sed -i -e "s/abbr1h/${PRE}abbr1h/g" $file
    sed -i -e "s/abbr0s/${PRE}abbr0s/g" $file
    sed -i -e "s/abbr1s/${PRE}abbr1s/g" $file
    sed -i -e "s/abbr1a/${PRE}abbr1a/g" $file    
    sed -i -e 's/#ifdef DEBUG/#ifdef DEBUGQ/g' $file
    # jede Zeile löschen, die #include "inline.h" enthält
    sed -i -e '/^#include "inline.h"/d' $file
done

echo "cd ${WORKINGDIR}/../${PROCDIR}/FormCalc_Virtuals/include"
cd ${WORKINGDIR}/../${PROCDIR}/FormCalc_Virtuals/include
PRE="${INITIAL1}_${FINAL}_"
for file in $(find . -type f -iname "${INITIAL1}_${FINAL}_vars.h"); do
    # alle Vorkommnisse von vars.h durch xxbar_vars.h ersetzen usw.
    sed -i -e "s/formfactors/${PRE}formfactors/g" $file
done
PRE="${INITIAL2}_${FINAL}_"
for file in $(find . -type f -iname "${INITIAL2}_${FINAL}_vars.h"); do
    # alle Vorkommnisse von vars.h durch xxbar_vars.h ersetzen usw.
    sed -i -e "s/formfactors/${PRE}formfactors/g" $file
done

# zurück an Anfang
cd ${WORKINGDIR}/

echo ""
echo "done."
echo "The only thing left to do now is to modify CalcRenConst.F and renconst.h ..."
echo "If you want to use finite Z/W-width, you should include this with "
echo "-cI*MW*WW in W-propagators or -cI*MZ*WZ in Z-propagators."
