# Epigenetic Regulation of Splicing
### Detecting Splicing Regulators via Correlation of RNA-seq and ChIP-seq ENCODE data
This program  uses publicly available data from ENCODE to detect epigenetic regulators of splicing of any human gene. This program correlates ChIP-seq peaks to RNA-seq transcripts at 500 bp bins 50 kb up and downstream of a designated gene, to assist in prediction of putative splicing regulators of circular and linear RNA. 

This program downloads all relevant ENCODE ChIP-seq and RNA-seq files, and outputs files that give the correlation between the peaks at ChIP-seq marks and RNA-seq quantification values across cell types (ENCODE Project Consortium 2012). There are three scripts needed to run this program. The first script (filter_metadata.py) only needs to be run once as it downloads an ENCODE metadata file, searches for relevant ChIP-seq and RNA-seq files for the analysis, and then downloads these files and produces a reduced metadata file. The second script (quant_and_bin.py) quantifies the RNA-seq files using salmon (Patro et al. 2017). and bins the ChIP-seq data with bedtools (Quinlan and Hall 2010), for a gene that the user specifies. This script produces a third script (do_analysis.sh) with all parameters preloaded to produce files that give the correlation between the ChIP-seq peaks and RNA-seq reads for a region along the genome. Therefore, the user only needs to run three scripts (Figure 1). This program is written using python, bash, and R. Required programs include bedtools and salmon; python modules include os, sys, stat, urllib, biomartpy (Smedley et al. 2015), subprocess, pandas, gzip, argparse, glob; and R packages include ggplot2 (Wickham 2009). 


## Steps
![userCommands](https://github.com/raungar/epi_reg/blob/master/readme_images/commands.png)

## Files and Folders
![folderStructure](https://github.com/raungar/epi_reg/blob/master/readme_images/action.png)
