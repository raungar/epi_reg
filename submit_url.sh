#!/bin/bash
#SBATCH --job-name="downloads_$file_name"
#SBATCH --output=jobs_output/downloaded.%j.out
#SBATCH --error=jobs_output/downloaded.%j.err
#SBATCH --time=45:00
#SBATCH --qos=normal
#SBATCH -p horence
#SBATCH --nodes=1
#SBATCH --mem=5Gb


python scripts/submit_url.py $1 $2
