#!/bin/bash

#SBATCH --job-name="analyze"
#SBATCH --output="jobs_output/analyze.%j.out"
#SBATCH --error="jobs_output/analyze.%j.err"
#SBATCH --time=30:00
#SBATCH --qos=normal
#SBATCH -p horence
#SBATCH --nodes=1
#SBATCH --mem=10Gb

#wrapper script to submit analysis

#give permissions to perform analysis
#needed because automatic generation does not give execution permissions
chmod u+x scripts/do_analysis.sh
./scripts/do_analysis.sh
