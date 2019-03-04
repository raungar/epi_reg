# Epigenetic Regulation of Splicing

### Step 1: Download all relevant files
##### filter_metadata.py
filter_metadata.py will create output directory with two subdirectories for the chip-seq and rna-seq downloads. This file will download the most recent metadata.tsv file from ENCODE and add the downloads to the appropriate directory. The metadata file of the downloaded files will be in the output file as reduced_metadata.tsv. Each rerun of this file should add any new files, though it will not update the existing files.


### Step 2: Quantification and Binning
##### salmon_compile_files.py
###### run_epireg.sh, rna_quant.py, chip_peaks.py
salmon_compile_files.py will perform quantification using salmon and find chip-seq peaks using bedtools. Any hg19 data will be converted to GRCh38. This file can work process single or paired end data. This is accomplished by distinguishing between rna/chip files and separately calling them as jobs (to speed up process) by calling run_epireg.sh. This calls either rna_quant.py or chip_peaks.py depending on the data to complete analysis.



