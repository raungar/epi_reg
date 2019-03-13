#!/bin/bash

#SBATCH --job-name="rerun-rna"
#SBATCH --output="jobs_output/rerun_rna.%j%.out"
#SBATCH --error="jobs_output/rerun_rna.%j%.err"
#SBATCH --time=1:15:00
#SBATCH --qos=normal
#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --mem=40Gb


echo $1

#export PATH=$PATH:~/local/numexpr/2.0.1/bin
#export PYTHONPATH=$PYTHONPATH:~/miniconda3/lib/python3.6/site-packages
#PYTHONPATH=$HOME/lib/python
#export PYTHONPATH=$PYTHONPATH:~/miniconda3/pkgs/biopython-1.72-py37h04863e7_0/lib/python3.7/site-packages/

if [[ "$1" == *"ChIP-seq"* ]]
then
	python scripts/chip_peaks.py $1 $2 $3 $4 $5 $6 $7
fi

if [[ "$1" == *"RNA-seq"* ]]
then
	python scripts/rna_quant.py $1 $2 $3 $4 $5 $6
fi

