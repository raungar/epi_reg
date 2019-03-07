#!/bin/bash

mkdir -p "output/quants/averaged/"

cell_types_list=`ls output/quants/RNA* | awk -F"_" '{print $3}' | sort | uniq`

for cell_type in $cell_types_list
do
	files=`ls output/quants/*$cell_type*`
	num_replicates=`ls output/quants/*$cell_type* | wc -l`
	if [[ $num_replicates == 1 ]]
	then
		file_name=`echo $files | awk -F"/" '{print $NF}'`
		cp $files output/quants/averaged/RNA_$cell_type.txt
	else
		declare -A transcripts

		for file in $files
		do
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

			done<$file
		done
		
		echo -e "ENST_ID\tSUM"> "output/quants/averaged/RNA_$cell_type.txt"

		for this_enst in "${!transcripts[@]}"
		do
			tpm_sum=${transcripts["$this_enst"]}
			tpm_ave=`echo "NA" | awk -v tpm_sum=$tpm_sum -v denom=$num_replicates '{ave=tpm_sum/denom; print ave}'`
			echo -e "$this_enst\t$tpm_ave" >> "output/quants/averaged/RNA_$cell_type.txt"
		done




	fi
		

done
