#!/bin/python

import sys
import os
import wget
import re
import multiprocessing as mp
import urllib.request
from multiprocessing.dummy import Pool
import subprocess


#this scripts downloads a metadata file of all human ENCODE experiments and then
#checks to see if it  passes filtering a filtering process. if so, it downloads relevant files
#in parallel by submitting many jobs. files are output in a rna/ directory for rna-seq files
#or in a chip/ directory for chip-seq files. reduced metadata file of all downloaded files is
#in reduced_metadata.tsv which will be used downstream

#to download files in parallel
def get_files(url):
	try:
		#record if it is chip-seq or rna-seq data in type
		if(re.search("bed",url)):
			type="chip"
		else:
			type="rna"
		#record output file name (type and id)
		file_name=str(type+"/"+url.split("/")[-1])
	
		#only download new files
		if not os.path.isfile(file_name):
			#submit download script via submit_url.sh
			os.system(str("sbatch scripts/submit_url.sh "+url+" "+file_name))
			log_file.write("Downloading: "+str(file_name)+"\n")

	except Exception as e:
		log_file.write("Download failed: " + str(url) + "\n")
		return 0

#downlaod files in parallel
def download(dic,type):
	url_list=[item[42] for item in dic.values()]
	pool = mp.Pool(processes=8)
	res = pool.map(get_files, url_list)
	log_file.write(str("Download complete: "+type+"\n"))

#check if line contains proper filtering requirements for chip-seq data
#returns 0 if it does not pass filtering requirements, zero if it does
def check_chip(line):
	passes=0 #set default to return 0 for not passing
	line_split=line.split("\t")

	#check if it is a bed narrowPeak file
	if(line_split[1] == "bed narrowPeak"):
		#check if the output type is "peaks"
		if(line_split[2]=="peaks"):
			#if matches h3k4me1, h3k4me3, h3k9me9, h3k27ac, k3k36me3 then it passes
			if(re.match("H3K(([27]{2}|[36]{2}|4|9)me3|4me1|27ac)",line_split[18])):
				passes=1
	return(passes)

#check if line contains proper filtering requirements for rna-seq data
#returns 0 if it does not pass filtering requirements, zero if it does
def check_rna(line):
	passes=0
	line_split=line.split("\t")
	#nested if statements to save time checking
	#check  if file format fastq
	if(line_split[1] == "fastq"):
		#check if output type is reads
		if(line_split[2]=="reads"):
			#check if biosample treatment is null
			if not (line_split[9]):
				#check if library depleted in rrna
				if(line_split[20]=="rRNA"):
					passes=1
	return(passes)

#read through the metadata file
#return everything passes the chip_seq thresholds as the first argument	
#return everything passes the rna_seq thresholds as the second argument	
def read_metadata(file_metadata):
	#make dictionaries
	rna_dic={}
	chip_dic={}

	try:
		#read file and determine if each line passes
		with open(file_metadata) as md_file:
			for line in md_file:
				line_split=line.split("\t")
			
				
				#if(line_split[50] != "" or line_split[51] != ""):
				#	continue

				#if audit errors, noncompliance, or action required skip
				if(line_split[51] != ""):
					continue
				#audit error (since it's the last line in the file it contains \n
				#so just search for any characters involved
				if(re.search("[a-zA-Z]",line_split[52])):
					continue

				#remove archived
				if(line_split[47] == "archived"):
					continue

				#get type
				if (line_split[4] == "ChIP-seq"):
					#if line passes chip-seq filters add to dic
					if(check_chip(line)==1):
						chip_dic[line_split[0]]=line_split
				elif (line_split[4] == "RNA-seq"):
					if(check_rna(line)==1):
						rna_dic[line_split[0]]=line_split
				else:
					continue
	except IOError:
		log_file.write("Error: Metadata File Not Found"+"\n")
	return(chip_dic,rna_dic)

#get which cells lines are the same in both the rna-seq and chip-seq data
#and return the dictionaries of only the files that have these cell line tyeps
def dics_same_cell_lines(chip_dic,rna_dic):
	#get the uniq chip and rna cell lines
	uniq_chip_cell_lines=set([col[6] for col in list(chip_dic.values())])
	uniq_rna_cell_lines=set([col[6] for col in list(rna_dic.values())])
	
	#get which cell lines were in both the uniq chip & rna cell lines
	same_cell_lines=uniq_chip_cell_lines.intersection(uniq_rna_cell_lines)
	chip_dic_final={}
	rna_dic_final={}

	#make a subsetted dictionary of only cell lines that both rna&chip contain
	for chip_key, chip_val in chip_dic.items():
		if chip_val[6] in same_cell_lines:
			chip_dic_final[chip_key]=chip_val
	for rna_key, rna_val in rna_dic.items():
		if rna_val[6] in same_cell_lines:
			rna_dic_final[rna_key]=rna_val
	return(chip_dic_final,rna_dic_final)

#downloads files and makes folders necessary for successful run
def make_paths_files():
	log_file.write("LOGFILE OUTPUT"+"\n")
	if not os.path.exists("chip"):
            os.mkdir("chip")
	if not os.path.exists("rna"):
            os.mkdir("rna")
	if not os.path.exists("metadata.tsv"):
		log_file.write("Downloading metadata.tsv..." +"\n")
		metadata_url=str("\"https://www.encodeproject.org/metadata/type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens/metadata.tsv\"")
		os.system(str("wget "+metadata_url)) #only downloads headers otherwise??
		while not os.path.exists("metadata.tsv"):
			time.sleep(0.2)	
			sys.stdout.write('.')
		log_file.write(" ---- Download complete"+"\n")


def main():
	
	#download files and make essential folders
	make_paths_files()
	file_metadata="metadata.tsv"

	#get dictionary of chip-seq and rna-seq files that pass filters
	chip_dic_all,rna_dic_all=read_metadata(file_metadata)
	#get dictionary of chip-seq and rna-seq files that contain overlapping cell lines
	chip_dic_reduced, rna_dic_reduced = dics_same_cell_lines(chip_dic_all,rna_dic_all)

	#download all resulting files in parallel 
	#will create many jobs
	download(chip_dic_reduced,"chip")
	download(rna_dic_reduced,"rna")
	
	#create a new file, reduced_metadata.tsv, that contains the metadata of the downloaded files
	reduced_metadata=open("reduced_metadata.tsv",'w')
	#print header
	with open('metadata.tsv') as md_file:
		header = md_file.readline()
		reduced_metadata.write(header)
	for reduced_vals_chip in chip_dic_reduced.values():
		reduced_metadata.write('\t'.join(map(str,reduced_vals_chip)))
	for reduced_vals_rna in rna_dic_reduced.values():
                reduced_metadata.write('\t'.join(map(str,reduced_vals_rna)))
	reduced_metadata.close()


log_file=open("log.txt",'w')
main()
log_file.close()

