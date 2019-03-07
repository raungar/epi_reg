#!/bin/bash

#make dir if does not exist
mkdir -p "output/matrix"


#declare associative arrays
declare -A marks
declare -A col_headers


#loop through output files to combine the ones that have the same mark
for f in `ls output/chip_peaks/ChIP*`
do 

	mark=`echo $f | awk -F"_" '{print $5}'`

	#if this is not the first instance of the mark add one to the number of times it's been found
	#else add it to the marks array and declare it to be found only once
	if test "${marks["$mark"]+isset}"
	then
		marks["$mark"]=$((${marks["$mark"]}+1))
	else
		marks["$mark"]=1
	fi

	
	#j is the number of times this mark has been found (first time is 1)
	j=${marks["$mark"]}
	prev_j=$(($j-1))


	#if it's the first time a mark has been found
	if [[ $j -eq 1 ]]
	then
		#define col_header to be mark and cell line
		col_header=`echo "$f" | awk -F"_" '{print $4"_"$5}'`
		#paste together the bed file such that the chrA loc1 loc2 --> chrA_loc1_loc2 is the first column
		#with the peak value from the mark file
		paste <(awk '{print $1"_"$2"_"$3}' output/bin_file_chr6_108509823_108734774.bed) <(awk '{print $4}' $f) > "output/matrix/temp_"$mark"_"$j".txt"
		#add to the col headers
		col_headers["$mark"]="location\t"$col_header
	else
		#define col header
		col_header=`echo "$f" | awk -F"_" '{print $4"_"$5}'`
		#if this mark has been found before, paste the column from this file to the previous file and save it in a new
		#temp file 
		paste "output/matrix/temp_"$mark"_"$prev_j".txt" <(awk '{print $4}' $f) > "output/matrix/temp_"$mark"_"$j".txt" 
		#add this column header to the array
		col_headers["$mark"]=${col_headers["$mark"]}"\t"$col_header
		#remove the previous temporary file since it's irrelevant
		rm "output/matrix/temp_"$mark"_"$prev_j".txt"

	fi
done
#head temp$j.txt


#for each mark found add a header column and create a final output matrix
for this_mark in "${!marks[@]}"
do
	#get teh file name and column header for this mark
	j=${marks[$this_mark]}
	file="temp_"$this_mark"_"$j".txt"
	this_col_header=${col_headers[$this_mark]}


	#actually add the header to the top of the file and save as a final file, remove the temporary file	
	echo -e "$this_col_header\n$(cat output/matrix/$file)" > output/matrix/chip_matrix_$this_mark.txt
	rm "output/matrix/temp_"$this_mark"_"$j".txt"

done
