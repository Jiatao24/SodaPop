#!/bin/bash

##### PARSE INPUT FROM C++ MAIN OR CMD LINE
OUT=$1
MAXGEN=$2
POPSIZE=$3
DT=$4
LONG=$5

##### CREATE DIRECTORIES FOR RESULTS
rm -r out/$OUT/barcodes; mkdir out/$OUT/barcodes
rm -r out/$OUT/graph; mkdir out/$OUT/graph
PREFIX=out/$OUT
echo $PREFIX

echo Extracting barcodes from $PREFIX/snapshots/...

#### CONVERT BINARY SNAPSHOTS TO TEXT FILES
FILES=$PREFIX/snapshots/*.snap
for filename in $FILES
do
	y=${filename%.001}
	sodasnap $filename $y.txt $LONG
done

#### EXTRACT AND SORT BARCODES
rm $PREFIX/avg_fitness.txt
FILES=$PREFIX/snapshots/*.snap.txt
for filename in $FILES
do
	y=${filename%%.txt}
	grep -w "C" $filename | cut -f1 | sort > $PREFIX/barcodes/${y##*/}.barcodes
	#### SUM POPULATION FITNESS FOR EACH TIME POINT AND DIVIDE BY POP SIZE
	grep -w "C" $filename | awk -v N=$3 '{sum += $4} END {print sum/N}' >> $PREFIX/avg_fitness.txt
done

echo Parsing unique barcodes...

#### PARSE UNIQUE BARCODES
FILES=$PREFIX/barcodes/*.barcodes
for filename in $FILES
do
	y=${filename%%.barcodes}
	uniq -c $filename | sed "s/^[ \t]*//" | awk -F' ' '{t = $1; $1 = $2; $2 = t; print; }' > $PREFIX/barcodes/${y##*/}.unique

	#### BREAK AT FIXATION POINT IF IT OCCURS
	if ! $(read -r && read -r)
	then
		#### OUTPUT FIXATION GENERATION TO FILE
	  	echo Fixation at ${y##*/} > $PREFIX/fixation.txt
	  	echo Fixation at ${y##*/}
	  	break
	fi < $PREFIX/barcodes/${y##*/}.unique
done

cat $PREFIX/barcodes/$OUT.gen0000000000.snap.unique > $PREFIX/barcodes/start.txt

echo Combining time series...

i=0
j=1

cat $PREFIX/barcodes/start.txt > $PREFIX/barcodes/series$i.txt

#### JOIN TIME FRAMES IN A SUITABLE FORMAT
# The format is for each line: BARCODE Gen(DTx0)_population Gen(DTx1)_population ...
for filename in `find $PREFIX/barcodes/ -maxdepth 1 -name '*.unique' | sort`
do
    # join:
    # -t' ' : ' ' is input/output separator
    # -e 0  : '0' replaces missing input fields
    # -a 1  : print unpairable lines from first file
    # -1 1  : join on field 1 of file 1 (this is default)
    # -2 1  : join on field 1 of file 2 (this is default)
    # -o 2.2 : output line format (file 2 field 2)

    # paste: merge lines of files
    # -d' ' : ' ' delimiter
    # Contents of previous series goes + new column of populations
    join -t' ' -e 0 -a 1 -1 1 -2 1 -o 2.2 $PREFIX/barcodes/series$i.txt $filename | paste -d' ' $PREFIX/barcodes/series$i.txt - > $PREFIX/barcodes/series$j.txt
    rm $PREFIX/barcodes/series$i.txt
    ((i++))
    ((j++))
done

# cut takes fields 1 and 3 to the end of the line
# Basically, snapshot from generation 0 is ignored.
cat $PREFIX/barcodes/series$i.txt | cut -d " " -f 1,3- > $PREFIX/ALL_generations.txt

rm $PREFIX/barcodes/series*.txt

#### PLOT RESULTS IN R SCRIPTS
# Presume that polyclonal_structure.R is in the same directory this script is in
location="$(dirname $(readlink -f ${BASH_SOURCE[0]}))"
Rscript "${location}/polyclonal_structure.R" /out/$OUT/ $DT

echo Done.
