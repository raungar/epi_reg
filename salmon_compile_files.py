#!/bin/python

import os
import sys
import stat
import urllib.request
#from pybiomart import Dataset
from biomart import BiomartServer
import biomartpy as bm
import pandas as pd
from pyliftover import LiftOver
import gzip


def make_bins(chr,start_init, end_init):
	start=str(start_init-50000)
	end=str(end_init+50000)
	
	os.system(str("echo "+chr)+"\t"+start+"\t"+end+str(" | sed 's/\s/\t/g' > output/chr_file_"+chr+"_"+start+"_"+end+".txt"))
	os.system(str("bedtools makewindows -b output/chr_file_"+chr+"_"+start+"_"+end+'''.txt -w 500 | awk '{fix=$2+1; print $1"\t"fix"\t"$3e}' > output/bin_file_'''+chr+"_"+start+"_"+end+".bed"))
	os.remove("output/chr_file_"+chr+"_"+start+"_"+end+".txt")


def make_folders():
	if not os.path.exists("output"):
		os.mkdir("output")

	if not os.path.exists("output/chip_peaks"):
		os.mkdir("output/chip_peaks")


	if not os.path.exists("output/hg19ToHg38.over.chain.gz"):
		urllib.request.urlretrieve("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", "output/hg19ToHg38.over.chain.gz")
		#urllib.request.urlretrieve("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver","output/liftOver")
		os.chmod("output/liftOver",777) #give permissions to operate on


	#MAYBE MOVE THE PATH AND FILE CHECKING TO MAIN GIRL OK 
	#(COULD MAYBE BE MESSED UP WHEN PARALLEL IMPLEMENTED)
	#build index if it does not exist
	if not os.path.exists("output/hs.grch38.index"):
		#download grch38 if file does not exist
		if not os.path.exists("output/Homo_sapiens.GRCh38.cdna.all.fa.gz"): 
			print("Downloading transcriptome fasta")
			os.system("wget -P output ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz")
			print("Download complete")
		#build salmon index
		print("building salmon index")
		os.system("salmon index -t output/Homo_sapiens.GRCh38.cdna.all.fa.gz -i output/hs.grch38.index")
		print("salmon index built")

	#create quants directory for salmon quantification if it odes not exist
	if not os.path.exists("output/quants"):
		os.mkdir("output/quants")




def main():
	print("HELLLO PLZ. PLZ :| ")
	cell_type_dic={}

	chr="chr6"
	start=108559823
	end=108684774
	gene="ENSG00000118689"

	make_folders()
	make_bins(chr,int(start),int(end))

	md_file=sys.argv[1]
	with open(md_file) as metadata:
		for md_line in metadata:
			md_line_split=md_line.split("\t")
			id=md_line_split[0]
			assay=md_line_split[4]
			cell_type=md_line_split[6]
			histone_mark=md_line_split[12].split("-")[0]
			paired_end=md_line_split[29]			
			pair2=md_line_split[30]
			assembly=md_line_split[37]
			
			if(cell_type in cell_type_dic.keys()):
				cell_type_dic[cell_type]+=1
			else:
				cell_type_dic[cell_type]=1


			if not(paired_end):
				R1=id
				R2=-1
			else:
				if(paired_end==1):
					R1=id
					R2=pair2
				else:
					continue

			rename_cell_type=cell_type.replace(" ", "-")
			rename_cell_type=rename_cell_type.replace("'","")
			print("THIS IS MY ASSAY: "+assay)
			if(assay=="ChIP-seq"):
				#OUTFILE, bed_id, chr, start, end, assemble
				outfile_name=str(assay+"_"+id+"_"+rename_cell_type+"_"+histone_mark+"_"+str(cell_type_dic[cell_type]))
				os.system("sbatch ./scripts/run_epireg.sh "+outfile_name+" "+id+" "+chr+" "+str(start)+" "+str(end)+" "+assembly)
				#os.system("sh scripts/run_epireg.sh "+id+" "+assay+" "+cell_type+" "+str(cell_type_dic[cell_type])+" "+histone_mark+" "+assembly)
				#run chip

			if(assay=="RNA-seq"):
				outfile_name=str(assay+"_"+id+"_"+rename_cell_type+"_"+str(cell_type_dic[cell_type]))
				##outfile, ensg_id, R1, R2
				os.system("sbatch ./scripts/run_epireg.sh "+outfile_name+" "+gene+" "+R1+" "+str(R2))
				#os.system("sh scripts/run_epireg.sh "+id+" "+assay+" "+cell_type+" "+str(cell_type_dic[cell_type])+" "+R1+" "+str(R2))
				#os.system("sh scripts/run_epireg.sh "+id+" "+assay+" "+cell_type+" "+str(cell_type_dic[cell_type])+" "+R1+" "+str(R2))
				#get rna
				#os.system(sh scripts/run_epireg.sh id assay cell_type cell_type_dic[cell_type] R1 R2)


	#	print(cell_type_dic)
			#print("\t".join([id,assay,cell_type,histone_mark,paired_end]))

	#OUTLINE: 
	#cell type dic
	#for line in metadata.tsv
	#get cell type
		#if file dont exist, plz create, name is just cell type
			#add cell type to dic, val is one
		#else, dic_num++ (counts # of replicates)
	#is rna
		#enst_id_dic=ensg_to_enst("ENSG00000118689")
		###get_rna_quant(list(enst_id_dic))

	#is chip
		#is hg19? hg19togrch38 that
		#is paired?
		#salmon
	#matchy matchy
		#below are col names k
		#bin deets (chr a:1,2) and quants from file
			#for each line, check if rna quants in there
				#if not, add zero for rna
				#if so, (FIGURE OUT WHEN MULTIPLE WHAT TO DO) print htat val
				#print dic val
				#print histone mark
				#print cell type
				
				#(WELLL RLLY>>>>> put in a line then print out the whole line u know)
		
	#for paired end, will need to only take the ones with a 1.
	#if can't find paired file, output an error message!
	#do_salmon(paired,"ENCFF091UZU","ENCFF163DLM")
	#do_salmon(paired,"ENCFF000HBA", "ENCFF000HBI")
#	hg19_to_GRCh38("ENCFF327LZT")

	#bed_file="ENCFF998UVS" #ENCFF327LZT
#	bed_file="ENCFF327LZT"
#	make_bins(chr,int(start), int(end)) #do before loop, should only need ot be done once
#	get_peaks(chr, int(start),int(end),bed_file)

main()
