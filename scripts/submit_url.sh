#!/bin/bash
#SBATCH --job-name="downloading_url"
#SBATCH --output=jobs_output/downloaded_url.%j.out
#SBATCH --error=jobs_output/downloaded_url.%j.err
#SBATCH --time=45:00
#SBATCH --qos=normal
#SBATCH -p horence
#SBATCH --nodes=1
#SBATCH --mem=10Gb

#submit script submit_url.py
python scripts/submit_url.py $1 $2

