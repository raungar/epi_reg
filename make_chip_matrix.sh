#!/bin/bash

outfolder=$1
chr=$2
start=$(($3-50000))
stop=$(($4+50000))
bed_file=$outfolder"/bin_file_"$chr"_"$start"_"$stop".bed"

awk '{print $1"_"$2"_"$3}' "$bed_file"

#make dir if does not exist
mkdir -p "$outfolder/matrix"

echo "YO"
echo $bed_file
echo "YOOO"


#declare associative arrays
declare -A marks
declare -A col_headers


#loop through output files to combine the ones that have the same mark
for f in `ls $outfolder"/chip_peaks/ChIP"*`
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
		col_header=`echo "$f" | awk -F"_" '{print $4}'`
		#paste together the bed file such that the chrA loc1 loc2 --> chrA_loc1_loc2 is the first column
		#with the peak value from the mark file
		paste <(awk '{print $1"_"$2"_"$3}' $bed_file ) <(awk '{print $4}' $f) > "$outfolder/matrix/temp_"$mark"_"$j".txt"
		#add to the col headers
		col_headers["$mark"]="location\t"$col_header
	else
		#define col header
		col_header=`echo "$f" | awk -F"_" '{print $4}'`
		#if this mark has been found before, paste the column from this file to the previous file and save it in a new
		#temp file 
		paste "$outfolder/matrix/temp_"$mark"_"$prev_j".txt" <(awk '{print $4}' $f) > "$outfolder/matrix/temp_"$mark"_"$j".txt" 
		#add this column header to the array
		col_headers["$mark"]=${col_headers["$mark"]}"\t"$col_header
		#remove the previous temporary file since it's irrelevant
		rm "$outfolder/matrix/temp_"$mark"_"$prev_j".txt"

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
	echo -e "$this_col_header\n$(cat $outfolder/matrix/$file)" > $outfolder/matrix/chip_matrix_$this_mark.txt
	rm "$outfolder/matrix/temp_"$this_mark"_"$j".txt"

done
#echo $marks

#uniq_marks=`printf '%s\n' $marks | awk -v RS='[[:space:]]+' '!a[$0]++{printf "%s%s", $0, RT}'`

#for this_mark in $uniq_marks 
#do
#	echo $this_mark
#	echo -e "$col_headers$this_mark\n$(cat $outfolder/matrix/temp_$this_mark_$j.txt)" > $outfolder/matrix/chip_matrix_$this_mark.txt
#done
#mv $outfolder/temp$j.txt $outfolder/matrix/chip_matrix.txt
