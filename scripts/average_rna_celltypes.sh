#!/bin/bash

outfolder=$1

#if does not exist, create dir averaged
mkdir -p "$outfolder/quants/averaged/"

#get list of all uniq cell types
cell_types_list=`ls $outfolder/quants/RNA* | awk -F"_" '{print $3}' | sort | uniq`

#create a file for each cell type that has the averaged value for that cell type for replicate cell types
#one line per enst where col1 is the enst name and col2 is the averaged value
for cell_type in $cell_types_list
do
	declare -A transcripts #make associative array of transcripts where key is enst and val is tpm value
	
	#get the list of files with this cell type
	#and record the number of replicates
	file_found=`ls $outfolder/quants/ |  grep -v "salmon" | grep $cell_type`
	num_replicates=`ls $outfolder/quants/ |  grep -v "salmon" | grep $cell_type | wc -l`
	files=$file_found

	#if it's the only file, then just take the enst name and val and save to quants/averaged
	if [[ $num_replicates == 1 ]]
	then
		file_name=`echo $files | awk -F"/" '{print $NF}'`
		#echo -e "ENST_ID\tSUM" > "$outfolder/quants/averaged/$cell_type.txt"
		cat "$outfolder/quants/$files" | awk '{print $1"\t"$4}' >  $outfolder/quants/averaged/$cell_type.txt
	#otherwise average the replicates
	else
		for x in "${!transcripts[@]}"; do printf "[%s]=%s\n" "$x" "${transcripts[$x]}" ; done

		#loop through all the replicate cell_type files
		for this_file in $files
		do
			this_file="$outfolder/quants/$this_file"
			while read line
			do

				enst=`echo $line | awk '{print $1}'` #get enst val
				tpm=`echo $line | awk '{print $4}'` #get tpm
				#see if this enst is already here, if so add tpm to the correct key
				if test "${transcripts["$enst"]+isset}"
				then
					val=`echo "NA" | awk -v sum="${transcripts["$enst"]}" -v tpm=$tpm '{new_sum=sum+tpm; print new_sum}'`
					transcripts["$enst"]="$val"
				#otherwise, just creat a new key/val pair
				else
					transcripts["$enst"]="$tpm"

				fi		

			done<"$this_file"
		done
		
		#once everything has been added, go back and divide by the number of replicates to actually get the average
		#and then print the result to quants/averaged/$cell_type.txt
		loop_i=0
		for this_enst in "${!transcripts[@]}"
		do
			loop_i=`echo $loop_i + 1 | bc`
			tpm_sum=${transcripts["$this_enst"]}
			tpm_ave=`echo "NA" | awk -v tpm_sum="$tpm_sum" -v denom="$num_replicates" '{ave=tpm_sum/denom; print ave}'`

			#overwrite if first time (allows for reruns)
			if [ $loop_i -eq 1 ]
			then
				echo -e "$this_enst\t$tpm_ave" > "$outfolder/quants/averaged/$cell_type.txt"
			else
				echo -e "$this_enst\t$tpm_ave" >> "$outfolder/quants/averaged/$cell_type.txt"
			fi
		done




	fi
		
	unset transcripts

done


#quants/enst: each file is an enst, one row per cell type where col1=cell type, col2=quant val
mkdir -p "$outfolder/quants/enst/"

line_i=0 #keep track of whether to write or overwrite (if first run, allows for rerun)

#get all averaged files
for f in `ls $outfolder/quants/averaged/*`
do
	#get cell types
	cell_type=`echo $f | awk -F"/" '{print $4}' | awk -F"." '{print $1}'`

	while read line
	do
		line_i=`echo $line_i + 1 | bc`
		enst=`echo $line | awk '{print $1}'` #get enst
		tpm=`echo $line | awk '{print $2}'` #get tpm
		#write the cell_type and tpm for each tpm
		if [[ $line_i -eq 1 ]]
		then
			echo -e "$cell_type\t$tpm" > "$outfolder/quants/enst/$enst.txt"

		else
			echo -e "$cell_type\t$tpm" >> "$outfolder/quants/enst/$enst.txt"
		fi
		 
	done<$f
done
