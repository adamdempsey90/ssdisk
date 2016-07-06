#!/bin/bash

DIR1=${1%/}
DIR2=${2%/}
NUM=$3
NEWNUM=${4:-$3}

cp $DIR1/gasdens$NUM.dat $DIR2/gasdens$NEWNUM.dat
cp $DIR1/gasvy$NUM.dat $DIR2/gasvy$NEWNUM.dat
cp $DIR1/gasvx$NUM.dat $DIR2/gasvx$NEWNUM.dat
cp $DIR1/gasenergy$NUM.dat $DIR2/gasenergy$NEWNUM.dat


cp $DIR1/domain_*.dat $DIR2/
cp $DIR1/*.var $DIR2/
cp $DIR1/*.par $DIR2/
cp $DIR1/planet*.dat $DIR2/
cp $DIR1/orbit*.dat $DIR2/
cp $DIR1/*_2d.dat $DIR2/
