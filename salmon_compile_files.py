#!/bin/python

import os
import sys
import stat
import urllib.request
from biomart import BiomartServer
import biomartpy as bm
import pandas as pd
from pyliftover import LiftOver
import gzip
import time
import argparse
import glob

def make_bins(chr,start_init, end_init,outfolder):
	start=str(start_init-50000)
	end=str(end_init+50000)
	
	os.system(str("echo "+chr)+"\t"+start+"\t"+end+str(" | sed 's/\s/\t/g' > "+outfolder+"/chr_file_"+chr+"_"+start+"_"+end+".txt"))
	os.system(str("bedtools makewindows -b "+outfolder+"/chr_file_"+chr+"_"+start+"_"+end+'''.txt -w 500 | awk '{fix=$2+1; print $1"\t"fix"\t"$3e}' > ''' +outfolder+"/bin_file_"+chr+"_"+start+"_"+end+".bed"))
	os.remove(outfolder+"/chr_file_"+chr+"_"+start+"_"+end+".txt")


def make_folders(outfolder):
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


	#MAYBE MOVE THE PATH AND FILE CHECKING TO MAIN GIRL OK 
	#(COULD MAYBE BE MESSED UP WHEN PARALLEL IMPLEMENTED)
	#build index if it does not exist
	if not os.path.exists("hs.grch38.index"):
		#download grch38 if file does not exist
		if not os.path.exists("Homo_sapiens.GRCh38.cdna.all.fa.gz"): 
			print("Downloading transcriptome fasta")
			os.system("wget ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz")
			print("Download complete")
		#build salmon index
		print("building salmon index")
		os.system("salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i hs.grch38.index")
		print("salmon index built")



def get_args():
	parser = argparse.ArgumentParser(prog='EPI_REG', usage='%(prog)s [options]')
	parser.add_argument("--metadata",type=str,help="metadata file to run on",nargs=1,required=True)
	parser.add_argument("--chr",type=str,help="chr name formatted like chr1",nargs=1,required=True)
	parser.add_argument("--start",type=int,help="start location",nargs=1,required=True)
	parser.add_argument("--end",type=int,help="end location",nargs=1,required=True)
	parser.add_argument("--gene",type=str,help="ENSG gene name",nargs=1,required=True)
	parser.add_argument("--outfolder",type=str,help="output folder",nargs=1,required=True)
	args=parser.parse_args()

	return(args)


def main():

	args=get_args()

	md_file=str(args.metadata[0])
	chr=str(args.chr[0])
	start=int(args.start[0])
	end=int(args.end[0])
	gene=str(args.gene[0])
	outfolder=str(args.outfolder[0])

	if not outfolder.isalnum():
		print("ERROR: OUTFOLDER MUST CONTAIN ONLY CHARACTERS OR NUMBERS")
		sys.exit()

	make_folders(outfolder)
	make_bins(chr,int(start),int(end),outfolder)


	cell_type_dic={}	
	line_count=-1 #start at -1 to account for header
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


			if((cell_type, histone_mark) in cell_type_dic.keys()):
				cell_type_dic[(cell_type,histone_mark)]+=1
			else:
				cell_type_dic[(cell_type,histone_mark)]=1


			if not(paired_end):
				R1=id
				R2=-1
			else:
				if(str(paired_end)=="1"):
					R1=id
					R2=pair2
				else:
					continue



			rename_cell_type=cell_type.replace(" ", "-")
			rename_cell_type=rename_cell_type.replace("'","")

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



	#wait until everything has finished
#	while((len(os.listdir(str(outfolder+"/chip_peaks")))+len(glob.glob(str(outfolder+"/quants/RNA*"))))<line_count):
#		print(".",end='')
#		time.sleep(5)

#	while((len(os.listdir("rna"))+len(os.listdir("rna"))<line_count):
#		time.sleep(5)

	new_script=open("scripts/do_analysis.sh","w")
	

	make_mx=str("./scripts/make_chip_matrix.sh "+outfolder+ " "+chr+" "+str(start)+" "+str(end))
	ave_rna=str("./scripts/average_rna_celltypes.sh "+outfolder)
	correlate=str("Rscript scripts/correlation.R "+outfolder)
	new_script.write("#!/bin/bash"+"\n"+make_mx+"\n"+ave_rna+"\n"+correlate+"\n")


main()
