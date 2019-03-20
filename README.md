# Epigenetic Regulation of Splicing
### Detecting Splicing Regulators via Correlation of RNA-seq and ChIP-seq ENCODE data
This program  uses publicly available data from ENCODE to detect epigenetic regulators of splicing of any human gene. This program correlates ChIP-seq peaks to RNA-seq transcripts at 500 bp bins 50 kb up and downstream of a designated gene, to assist in prediction of putative splicing regulators of circular and linear RNA. 

This program downloads all relevant ENCODE ChIP-seq and RNA-seq files, and outputs files that give the correlation between the peaks at ChIP-seq marks and RNA-seq quantification values across cell types (ENCODE Project Consortium 2012). There are three scripts needed to run this program. The first script (filter_metadata.py) only needs to be run once as it downloads an ENCODE metadata file, searches for relevant ChIP-seq and RNA-seq files for the analysis, and then downloads these files and produces a reduced metadata file. The second script (quant_and_bin.py) quantifies the RNA-seq files using salmon (Patro et al. 2017). and bins the ChIP-seq data with bedtools (Quinlan and Hall 2010), for a gene that the user specifies. This script produces a third script (do_analysis.sh) with all parameters preloaded to produce files that give the correlation between the ChIP-seq peaks and RNA-seq reads for a region along the genome. Therefore, the user only needs to run three scripts (Figure 1). This program is written using python, bash, and R. Required programs include bedtools and salmon; python modules include os, sys, stat, urllib, biomartpy (Smedley et al. 2015), subprocess, pandas, gzip, argparse, glob; and R packages include ggplot2 (Wickham 2009). 


## Steps

#### 1. Download Files    
###### Create a new directory, and make scripts subdirectory with all the scripts inside it. This will download all relevant files in parallel if using SLURM. Each time this is rerun, it will add new files that have been added to ENCODE, but will not update pre-existing files.  
This script automatically downloads ENCODE human metadata from https://www.encodeproject.org/metadata/type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens/metadata.tsv which is saved to metadata.tsv . ChIP-seq files that are in “bed narrowPeak” file format with have an output type of peaks and histone marks that include H3K4me1, H3K4me3, H3K9me3, H3K27ac, H3K36me3 are selected. RNA-seq files that are in fastq file format with an output type of reads and no biosample treatment and that are depleted in rRNA are retained. Any files with archived data and files with messages about audit errors and non-compliance are discarded. Only experiments performed on cell types in both the RNA-seq and ChIP-seq files are kept. Once the final subset of files has been determined, a metadata file of only the final files is then created as reduced_metadata.tsv. The ChIP-seq files are downloaded into a folder named “chip” and RNA- seq files are downloaded into a folder named “rna” in a parallel fashion. This is done by submitting a job via submit_url.sh for each file download by calling the file submit_url.py within filter_metadata.py. Note that this submit many jobs. Each time the script is run, it will only download files if the files do not exist, and therefore if the jobs have timed-out or been limited, the user can simply rerun the program. Removing “metadata.tsv” and rerunning the script can act as a way to update metadata.tsv and add more files, as any additional files added will be added to the appropriate directories, but not that this will not update existing files.
    
#### 2. Quantification and Binning    
###### This performs RNA-seq quantification in parallel if it does not already exist, and bins the ChIP-seq peaks for the relevant genes. If this has already been run, use the nonparallel_quant_and_bin.py file instead.
*quant_and_bin.py* The script takes user input of the reduced metadata file, the gene name, the chromosome name, the start location, the end location, and the name of the output folder. All future files will be outputted into this output folder. An optional index can be provided, otherwise hg38 index is used. If known, the user can specify a preferred ENST transcript to include as an optional parameter.The script performs error checking to make sure that the user is inputting the correct parameters. This script then creates windows for the specified location using bedtools and outputs this file as bin_file_[chr]_[start]_[stop].bed. To save time and parallelize the process, the metadata file is read in a job is created to either perform the RNA- seq quantification or the ChIP-seq peak calling. This is done by submitting run_epireg.sh. This will script runs chip_peaks.py if the file is a ChIP-seq file or rna_quants.py if the file is a RNA-seq file.     
*chip_peaks.py* The file first checks if the assembly is hg19 or GRCh38. Currently hg19 data is not supported, though this will be implemented in future versions, so such files are ignored. If this passes, the bedtools is used to obtain the peaks only at the region within 50kb up and downstream of the start and stop respectively. This is outputted in outfolder/chip_peaks.     
*rna_quant.py* This file uses salmon to quantify RNA-seq reads. Salmon was used rather than kallisto because kallisto can only properly deal with single-end data if experimental fragment data is known. This information is unlisted in the metadata, and a kallisto developer suggested to use salmon as an alternative in this situation (Chapman). The script determines if the data is single-end or paired-end, and checks to see if a salmon quantification file already exists. If it does not exist salmon quantification is performed and outputted in the salmon directory. The script then queries biomart using the given ENSG name to retrieve proper ENST names, which is then used to get only relevant transcripts from the salmon quantification. All relevant transcripts and quantifications are stored with an informative filename of the cell type, replicate number, and sample id in output/quants directory.    


#### 3. Perform Analysis    
###### This gives summary data files of ChIP-seq data in outfolder/matrix and RNA-seq data in outfolder/quants/enst & outfolder/quants/averaged. This produces the ultimate correlation files in outfolder/correlation. This script is automatically generated after running step 2. 
*do_analysis.sh* Once all jobs have been completed, analyze.sh acts as a wrapper scripts then calls make_chip_matrix.sh, average_rna_celltypes.sh, and correlation.R. This script is automatically generated with correct parameters.
*make_chip_matrix.sh* This script makes one matrix per histone mark from the chip peak files and outputs it into outfolder/matrix. The first column is the bin location, the header is the cell type, and the matrix values are the peak at that location and cell-type.     
*average_rna_celltypes.sh* This script takes all the RNA quantification files and creates two directories that can be used for different interpretation of the data. The directory outfolder/quants/enst contains a file of ENST transcript and the first column is the cell type and the second column is the average transcript product at that location. The directory outfolder/quants/averaged contains one file per cell type that lists in the first column the ENST transcript and the second column contains the average transcript product at this location.     
*correlations.R*  This script is called to produce correlation files which are outputted in the outfolder/correlation directory. For all combinations of histone marks and transcripts, a correlation file is produced between with the columns being the bin location, pearson correlation, and spearman correlation. This averages any duplicate ChIP-seq cell types. This can be used to identify regions with high correlation between ChIP-seq peaks and a gene transcript.

![userCommands](https://github.com/raungar/epi_reg/blob/master/readme_images/commands.png)
    



## Directory Structure and Graphical Representation of Script
![folderStructure](https://github.com/raungar/epi_reg/blob/master/readme_images/action.png)
