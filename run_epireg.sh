#!/bin/bash

#SBATCH --job-name="epireg_$1"
#SBATCH --output=/scratch/PI/horence/rachel/epi_reg/job_output/epireg_$1.%j.out
#SBATCH --error=/scratch/PI/horence/rachel/epi_reg/job_output/epireg_$1.%j.err
#SBATCH --time=30:00
#SBATCH --qos=normal
#SBATCH -p horence
#SBATCH --nodes=1
#SBATCH --mem=40Gb

if [[ "$1" == *"ChIP-seq"* ]]
then
	chip_peaks.py $1 $2 $3 $4 $5 $6
fi

if [[ "$1" == *"RNA-seq"* ]]
then
	python rna_quant.py $1 $2 $3 $4 
fi
