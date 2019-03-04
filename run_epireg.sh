#!/bin/bash

#SBATCH --job-name="epireg_$1"
#SBATCH --output=/scratch/PI/horence/rachel/epi_reg/job_output/epireg_$1.%j.out
#SBATCH --error=/scratch/PI/horence/rachel/epi_reg/job_output/epireg_$1.%j.err
#SBATCH --time=30:00
#SBATCH --qos=normal
#SBATCH -p horence
#SBATCH --nodes=1
#SBATCH --mem=40Gb

echo "THIS IS $1"

if [[ "$1" == *"ChIP-seq"* ]]
then
	#$1 $2 $3 $4 $5 $6 $7
	#OUTFILE, bed_id, chr, start, end, assemble
	echo "CHIP"
fi

if [[ "$1" == *"RNA-seq"* ]]
then
	#$1 $2 $3 $4 $5 $6 $7
	#outfile, ensg_id, R1, R2
	echo "RNA"
fi
