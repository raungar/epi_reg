#!/bin/bash

outfolder=$1

mkdir -p "$outfolder/quants/averaged/"

cell_types_list=`ls $outfolder/quants/RNA* | awk -F"_" '{print $3}' | sort | uniq`

for cell_type in $cell_types_list
do
	file_found=`ls $outfolder/quants/ |  grep -v "salmon" | grep $cell_type`
	num_replicates=`ls $outfolder/quants/ |  grep -v "salmon" | grep $cell_type | wc -l`
	files=$file_found


	if [[ $num_replicates == 1 ]]
	then
		file_name=`echo $files | awk -F"/" '{print $NF}'`
		#echo -e "ENST_ID\tSUM" > "$outfolder/quants/averaged/$cell_type.txt"
		cat "$outfolder/quants/$files" | awk '{print $1"\t"$4}' >  $outfolder/quants/averaged/$cell_type.txt
	else
		declare -A transcripts

		for this_file in $files
		do
			#echo "$file"
			this_file="$outfolder/quants/$this_file"
			#echo $this_file
			while read line
			do

				enst=`echo $line | awk '{print $1}'`
				tpm=`echo $line | awk '{print $4}'`
				if test "${transcripts["$enst"]+isset}"
				then
					val=`echo "NA" | awk -v sum="${transcripts["$enst"]}" -v tpm=$tpm '{new_sum=sum+tpm; print new_sum}'`
	
					#echo -e "in if: $enst\t$val"
					transcripts["$enst"]="$val"
				else
					transcripts["$enst"]="$tpm"
					#echo -e "in else: $enst\t$tpm"

				fi		

			done<"$this_file"
		done
		

		#echo -e "ENST_ID\tSUM" > "$outfolder/quants/averaged/$cell_type.txt"
		loop_i=0
		for this_enst in "${!transcripts[@]}"
		do
			loop_i=`echo $loop_i + 1 | bc`
			tpm_sum=${transcripts["$this_enst"]}
			tpm_ave=`echo "NA" | awk -v tpm_sum=$tpm_sum -v denom=$num_replicates '{ave=tpm_sum/denom; print ave}'`
			if [ $loop_i -eq 1 ]
			then
				echo -e "$this_enst\t$tpm_ave" > "$outfolder/quants/averaged/$cell_type.txt"
			else
				echo -e "$this_enst\t$tpm_ave" >> "$outfolder/quants/averaged/$cell_type.txt"
			fi
		done




	fi
		

done

mkdir -p "$outfolder/quants/enst/"

line_i=0
for f in `ls $outfolder/quants/averaged/*`
do
	cell_type=`echo $f | awk -F"/" '{print $4}' | awk -F"." '{print $1}'`

	#line_i=$(($line_i+1)
	while read line
	do
		line_i=`echo $line_i + 1 | bc`
		echo $line_i
		enst=`echo $line | awk '{print $1}'`
		tpm=`echo $line | awk '{print $2}'`
		if [[ $line_i -eq 1 ]]
		then
			echo -e "$cell_type\t$tpm" > "$outfolder/quants/enst/$enst.txt"

		else
			echo -e "$cell_type\t$tpm" >> "$outfolder/quants/enst/$enst.txt"
		fi
		 
	done<$f
done



