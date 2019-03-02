#!/bin/python

import os
import stat
import urllib.request
from biomart import BiomartServer
import biomartpy as bm
import pandas as pd

def do_salmon(paired,R1,R2):
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

	#perform salmon quantification on single end or paired end data
	print("salmon quanitification")
	if(paired=="False"):
		print("single end salmon")
		os.system(str("salmon quant -i output/hs.grch38.index -l A -r output/rna/"+R1+".fastq.gz -p 8 -o output/quants/salmon_"+R1))
	else:
		print("paired end salmon")
		os.system(str("salmon quant -i output/hs.grch38.index -l A -1 output/rna/"+R1+".fastq.gz -2 output/rna/"+R2+".fastq.gz -p 8 -o output/quants/salmon_"+R1+"_"+R2))


#make bins for chip peak analysis
def make_bins(chr,start_init, end_init):
	start=str(start_init-50000)
	end=str(end_init+50000)
	
	os.system(str("echo "+chr)+"\t"+start+"\t"+end+str(" | sed 's/\s/\t/g' > output/chr_file_"+chr+"_"+start+"_"+end+".txt"))
	os.system(str("bedtools makewindows -b output/chr_file_"+chr+"_"+start+"_"+end+'''.txt -w 500 | awk '{fix=$2+1; print $1"\t"fix"\t"$3e}' > output/bin_file_'''+chr+"_"+start+"_"+end+".bed"))
	os.remove("output/chr_file_"+chr+"_"+start+"_"+end+".txt")


#get peak score
def get_peaks(chr,start_init, end_init,chip_file_id):
	#region is 50,000kb +/- of start and stop
	start=str(start_init-50000)
	end=str(end_init+50000)

	#if the directory does not exist, create it
	if not os.path.exists("output/chip_peaks"):
		os.mkdir("output/chip_peaks")
	
	#get peaks at intersection with bin file
	os.system(str("bedtools intersect -wao -a  output/bin_file_"+chr+"_"+start+"_"+end+".bed -b output/rna/"+chip_file_id+".bed.gz > output/chip_peaks/"+chip_file_id+".bed"))
	#sort for merge step
	os.system(str("sort -k1,1 -k2,2n output/chip_peaks/"+chip_file_id+".bed > output/chip_peaks/sorted_"+chip_file_id+".bed"))
	os.system(str("bedtools merge -i output/chip_peaks/sorted_"+chip_file_id+".bed -c 14 -o max > output/chip_peaks/merged_"+chip_file_id+".bed"))
	os.remove(str("output/chip_peaks/sorted_"+chip_file_id+".bed"))
	os.remove(str("output/chip_peaks/"+chip_file_id+".bed"))



def hg19_to_GRCh38(file):
	if not os.path.exists("output/hg19ToHg38.over.chain.gz"):
		urllib.request.urlretrieve("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", "output/hg19ToHg38.over.chain.gz")
		urllib.request.urlretrieve("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver","output/liftOver")
		os.chmod("output/liftOver",777) #give permissions to operate on
	os.system(str("output/liftOver output/rna/"+file+".bed.gz output/hg19ToHg38.over.chain.gz output/rna/lifted_"+file+".bed.gz unMapped"))


#gets transcript ids and location for ensembl gene name
#and returns this information as a dictionary (keys=enst_id,vals=information list)
def ensg_to_enst(gene_id):
	#define biomart parameters
	mart_name = "ensembl"
	dataset = "hsapiens_gene_ensembl"
	attributes=["ensembl_transcript_id","chromosome_name","start_position","end_position"]
	filters={"ensembl_gene_id":[gene_id]}
	#using biomartpy, obtain pandas df of requested information
	biomart_df = bm.make_lookup(mart_name=mart_name, dataset=dataset, attributes=attributes, filters=filters)
	#transform this data frame to a dictionary where the row names (enst_id) are the keys
	#and the values are the rest of the requested information as a list
	enst_dic=biomart_df.to_dict("index")
	return(enst_dic)

def get_rna_quant(enst_id):
	#get number of max matches so once all ids are found they can stop running through them
	my_file="output/quants/salmon_ENCFF091UZU_ENCFF163DLM/quant.sf"
	with open(my_file, "r") as quant_file:
		for quant_line in quant_file:
			found=-1
			for this_id in enst_id:
				if(this_id) in quant_line:
					found=this_id
			#if any(this_id in quant_line for this_id in enst_id):
			if(found !=-1 ):
				print(found)
				enst_id.remove(found) #remove to save time
				print(quant_line.strip())
				#if all enst_id are found, stop searching
				if(len(enst_id)==0):
					break				
def combine_files(start,stop):
	rpkm=100*tpm/(stop-start)
	

def main():
	#os.system("conda activate salmon")

	#os.system("source ~/bisque/bqenv/bin/activate")

	enst_id_dic=ensg_to_enst("ENSG00000118689")

	print(list(enst_id_dic))

	get_rna_quant(list(enst_id_dic))
	#s=Estimated average fragment length
	#l=Estimated average fragment length
	#these are kallisto parameters, required for single end.
	#if paired end kallisto calculates this, pass -1, -1 to symbolize paired end
	s=200
	l=10
	#do_kallisto(s,l)
	paired="True"
	#do_salmon(paired, "ENCFF163DLM", -1)
	#do_salmon(paired, "ENCFF001RMC,-1)

	#for paired end, will need to only take the ones with a 1.
	#if can't find paired file, output an error message!
	#do_salmon(paired,"ENCFF091UZU","ENCFF163DLM")
	#do_salmon(paired,"ENCFF000HBA", "ENCFF000HBI")
	#hg19_to_GRCh38("ENCFF327LZT")

	chr="chr6"
	start=108559823
	end=108684774
	#bed_file="ENCFF998UVS" #ENCFF327LZT
	bed_file="ENCFF327LZT"
#	make_bins(chr,int(start), int(end)) #do before loop, should only need ot be done once
#	get_peaks(chr, int(start),int(end),bed_file)

main()
