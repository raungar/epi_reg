#!/bin/bash

#SBATCH --job-name="run_chip_or_rna"
#SBATCH --output="jobs_output/rerun_rna.%j%.out"
#SBATCH --error="jobs_output/rerun_rna.%j%.err"
#SBATCH --time=1:15:00
#SBATCH --qos=normal
#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --mem=40Gb



#check if it is a rna-seq or chip-seq file, and then call appropriate script
if [[ "$1" == *"ChIP-seq"* ]]
then
	python scripts/chip_peaks.py $1 $2 $3 $4 $5 $6 $7
fi

if [[ "$1" == *"RNA-seq"* ]]
then
	python scripts/rna_quant.py $1 $2 $3 $4 $5 $6
fi

