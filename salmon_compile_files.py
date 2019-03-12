#!/bin/python

import os
import sys
import stat
import urllib.request
from pyliftover import LiftOver
import gzip
import argparse
import glob

#make a bin file for bedtools that start 50kb up and downstream of start and stop
#with bin size of 500 bp. do here since only needs to be done once
def make_bins(chr,start_init, end_init,outfolder):
	#50kb up and downstream
	start=str(start_init-50000)
	end=str(end_init+50000)
	
	#create naming
	os.system(str("echo "+chr)+"\t"+start+"\t"+end+str(" | sed 's/\s/\t/g' > "+outfolder+"/chr_file_"+chr+"_"+start+"_"+end+".txt"))
	#uses bedtools makewindows to create location windows
	os.system(str("bedtools makewindows -b "+outfolder+"/chr_file_"+chr+"_"+start+"_"+end+'''.txt -w 500 | awk '{fix=$2+1; print $1"\t"fix"\t"$3e}' > ''' +outfolder+"/bin_file_"+chr+"_"+start+"_"+end+".bed"))
	os.remove(outfolder+"/chr_file_"+chr+"_"+start+"_"+end+".txt")

#make necessary folders if they do not exist
#and download necessary files if they do not exist
#and build salmon quantification file if it does not exist
def make_folders(outfolder, custom_index):
	if not os.path.exists(outfolder):
		os.mkdir(outfolder)
	if not os.path.exists("jobs_output"):
		os.mkdir("jobs_output")

	if not os.path.exists(str(outfolder+"/chip_peaks")):
		os.mkdir(str(outfolder+"/chip_peaks"))
	#create quants directory for salmon quantification if it odes not exist
	if not os.path.exists(str(outfolder+"/quants")):
		os.mkdir(str(outfolder+"/quants"))


	if not os.path.exists("/hg19ToHg38.over.chain.gz"):
		urllib.request.urlretrieve("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", "hg19ToHg38.over.chain.gz")
		#urllib.request.urlretrieve("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver","liftOver")
		#os.chmod("liftOver",0o777) #give permissions to operate on


	#build index if it does not exist
	if not os.path.exists("hs.grch38.index") and custom_index=-"-1":
		#download grch38 if file does not exist
		if not os.path.exists("Homo_sapiens.GRCh38.cdna.all.fa.gz"): 
			print("Downloading transcriptome fasta")
			os.system("wget ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz")
			print("Download complete")
		#build salmon index
		print("building salmon index")
		os.system("salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i hs.grch38.index")
		print("salmon index built")
	else:
		os.system("salmon index -t " +custom_index+ " -i hs.grch38.index")

#take user input, force all arguments (no optional), and return args
def get_args():
	parser = argparse.ArgumentParser(prog='EPI_REG', usage='%(prog)s [options]')
	parser.add_argument("--metadata",type=str,help="metadata file to run on",nargs=1,required=True)
	parser.add_argument("--chr",type=str,help="chr name formatted like chr1",nargs=1,required=True)
	parser.add_argument("--start",type=int,help="start location",nargs=1,required=True)
	parser.add_argument("--end",type=int,help="end location",nargs=1,required=True)
	parser.add_argument("--gene",type=str,help="ENSG gene name",nargs=1,required=True)
	parser.add_argument("--outfolder",type=str,help="output folder",nargs=1,required=True)
	parser.add_argument("--custom_index",type=str,help="full path to custom index",nargs=1,required=False)

	args=parser.parse_args()

	return(args)


def main():

	print("Beginning Analysis")
	args=get_args()

	#save args into associated variable name
	md_file=str(args.metadata[0])
	chr=str(args.chr[0])
	start=int(args.start[0])
	end=int(args.end[0])
	gene=str(args.gene[0])
	outfolder=str(args.outfolder[0])
	#check if custom index exists
	if not str(args.custom_index[0]):
		custom_index="-1"
	else:		
		custom_index=str(args.custom_index[0])


	#do not allow outfolder name to not be alphanumeric
	#becuase could potentially break downstream awk-ing
	if not outfolder.isalnum():
		print("ERROR: OUTFOLDER MUST CONTAIN ONLY CHARACTERS OR NUMBERS")
		sys.exit()

	#make all necessary folders and create all necessary files
	make_folders(outfolder,custom_index)

	#make the bin file that will be used by bedtools for binning
	#do here since only needs to be made once
	make_bins(chr,int(start),int(end),outfolder)


	cell_type_dic={} #count how many times this cell type has been seen (key=(cell_type,mark) and val=count)
	line_count=-1 #start at -1 to account for header
	#open metadata file and submit job per file
	with open(md_file) as metadata:
		for md_line in metadata:
			line_count+=1
			md_line_split=md_line.split("\t")
			id=md_line_split[0]
			assay=md_line_split[4]
			cell_type=md_line_split[6]
			histone_mark=md_line_split[18].split("-")[0]
			paired_end=md_line_split[35]			
			pair2=md_line_split[36]
			assembly=md_line_split[43]

			#print("HISTONE MARK: "+histone_mark + " , assembly: "+assembly+" , celltype: "+cell_type)
			#print("paired_end: "+paired_end+", pair2: "+pair2, "assay: "+assay)

			#record count of how many times this celltype/histone mark pair has been seen
			if((cell_type, histone_mark) in cell_type_dic.keys()):
				cell_type_dic[(cell_type,histone_mark)]+=1
			else:
				cell_type_dic[(cell_type,histone_mark)]=1


			#if the data is not paired end, R2 should be set to -1 to inform
			if not(paired_end):
				R1=id
				R2=-1
			else:
				#if it's the first in the pair do analysis
				#otherwise skip (to avoid doing repeat analysis)
				if(str(paired_end)=="1"):
					R1=id
					R2=pair2
				else:
					continue


			#fix cell type naming
			rename_cell_type=cell_type.replace(" ", "-")
			rename_cell_type=rename_cell_type.replace("'","")

			#call scripts/run_epireg with proper parameters
			#will submit a job that is either chip_peaks.py for chip-seq to bin files with bedtools
			#or rna_quant.py for rna-seq to perform quantification with salmon
			if(assay=="ChIP-seq"):
				#OUTFILE, bed_id, chr, start, end, assemble
				outfile_name=str(outfolder+"/chip_peaks/"+assay+"_"+id+"_"+rename_cell_type+"_"+histone_mark+"_"+str(cell_type_dic[(cell_type,histone_mark)]))
				os.system("sbatch ./scripts/run_epireg.sh "+outfile_name+" "+id+" "+chr+" "+str(start)+" "+str(end)+" "+assembly+" "+outfolder)
				#os.system("sh scripts/run_epireg.sh "+id+" "+assay+" "+cell_type+" "+str(cell_type_dic[cell_type])+" "+histone_mark+" "+assembly)
				#run chip

			if(assay=="RNA-seq"):
				outfile_name=str(outfolder+"/quants/"+assay+"_"+id+"_"+rename_cell_type+"_"+str(cell_type_dic[(cell_type,histone_mark)]))
				##outfile, ensg_id, R1, R2
				os.system("sbatch ./scripts/run_epireg.sh "+outfile_name+" "+gene+" "+R1+" "+str(R2)+" "+outfolder)
				#os.system("sh scripts/run_epireg.sh "+id+" "+assay+" "+cell_type+" "+str(cell_type_dic[cell_type])+" "+R1+" "+str(R2))
				#os.system("sh scripts/run_epireg.sh "+id+" "+assay+" "+cell_type+" "+str(cell_type_dic[cell_type])+" "+R1+" "+str(R2))
				#get rna
				#os.system(sh scripts/run_epireg.sh id assay cell_type cell_type_dic[cell_type] R1 R2)



	#write the script do_analysis.sh based on inputs here
	#this script calls three scripts to perform analyses sequentially
	new_script=open("scripts/do_analysis.sh","w")	

	#make_chip_matrix makes a matrix per histone mark where the cols are cell lines, the rows 
	#are bin locations, and the values are peak values
	make_mx=str("./scripts/make_chip_matrix.sh "+outfolder+ " "+chr+" "+str(start)+" "+str(end))

	#average rna_cell_types.sh produces two folders with different version of the same content
	#quants/averaged: each file is a cell type, one row per enst where col1=enst, col2=quant val
	#quants/enst: each file is an enst, one row per cell type where col1=cell type, col2=quant val
	ave_rna=str("./scripts/average_rna_celltypes.sh "+outfolder)

	#produces correlation file per histone mark x enst
	#col1=bin locations, col2=pearson, col3=spearman correlation
	correlate=str("Rscript scripts/correlation.R "+outfolder)

	new_script.write("#!/bin/bash"+"\n"+make_mx+"\n"+ave_rna+"\n"+correlate+"\n")


main()
