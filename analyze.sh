#!/bin/bash


outfolder=$1
if [ -z "$1" ]
then
	echo "No argument supplied. Please use an argument that is the output folder to analyze"
	exit 1
fi



if [ -d $outfolder ]
then
	./scripts/make_chip_matrix.sh "$outfolder"
	./scripts/average_rna_celltypes.sh "$outfolder"
	Rscript scripts/correlation.R "$outfolder"

else
	echo "Outfolder does not exist. Please use an argument that is the output folder to analyze"
	exit 1
fi
