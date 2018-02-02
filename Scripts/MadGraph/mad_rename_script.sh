#!/bin/bash
# Copyright (C) Matthias Kesenheimer - All Rights Reserved
# Written by Matthias Kesenheimer <m.kesenheimer@gmx.net>, 2017

# dieses Skript ersetzt alle in Madgraph vorkommenden Variablen MB1 und MB2
# durch MBL und MBR. Diese Umbenennung wurde vorgenommen um Konflikte mit
# anderen Variablen zu umgehen und um die Commonblöcke nicht ändern zu 
# müssen. Dieses Skript muss im gleichen Ordner der zu Verändernden Dateien
# ausgeführt werden.
# Dieses Skript ist mit Vorsicht zu genießen.
for x in $(find . -type f -iname "*.f" -o -iname "*.inc"); do
    sed -i -e 's/\<MB1\>/MBL/g' $x
    sed -i -e 's/\<MB2\>/MBR/g' $x
done
