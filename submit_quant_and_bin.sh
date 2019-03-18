#!/bin/bash
#SBATCH --job-name="hotair"
#SBATCH --output=jobs_output/quantbin_fat.%j.out
#SBATCH --error=jobs_output/quantbin_fat.%j.err
#SBATCH --time=1:30:00
#SBATCH --qos=normal
#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --mem=30Gb

#wrapper script to submit quant_and_bin.py

python scripts/quant_and_bin.py --metadata "reduced_metadata.tsv" --chr "chr12"  --start 53962308 --end 53974956   --gene "ENSG00000228630" --outfolder "hotair"

