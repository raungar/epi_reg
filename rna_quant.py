#!/bin/python

import os
import sys
import urllib
import biomartpy as bm

#this files performs salmon quantification if not already done for the input file
#then gets the relevant ensg --> enst and adds on any optional enst provided
#then gets the relevant transcripts from the quantification files and outputs them
#in the folder outfolder/quants/

#run salmon quantification
def do_salmon(R1,R2):
	#create quants directory for salmon quantification if it odes not exist
	if not os.path.exists("salmon"):
		os.mkdir("salmon")

	#perform salmon quantification on single end or paired end data
	if(str(R2)==str("-1")):
		print("single end salmon")
		file_name=str("salmon_"+R1)
		os.system(str("salmon quant -i hs.grch38.index -l A -r rna/"+R1+".fastq.gz -p 8 -o salmon/salmon_"+R1))
	else:
		file_name=str("salmon_"+R1+"_"+R2)
		print("paired end salmon")
		os.system(str("salmon quant -i hs.grch38.index -l A -1 rna/"+R1+".fastq.gz -2 rna/"+R2+".fastq.gz -p 8 -o salmon/salmon_"+R1+"_"+R2))
	return(file_name)

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


#for the salmon quantification, find the enst transcripts and print them to a 
#file with a more descriptive name (outfile) for future analysis
def get_rna_quant(enst_id,rna_file,outfolder,outfile):
	if not os.path.exists(str(outfolder+"/quants")):
		os.mkdir(str(outfolder+"/quants"))

	print(outfile)
	#salmon file to look through
	#my_file=str("salmon/"+rna_file+"/quant.sf")

	#open rna_output as file to write to
	with open(outfile,"w") as rna_output:
		#open rna output file
		with open(rna_file, "r") as quant_file:
			for quant_line in quant_file:
				found=-1 #keeps track of whether the enst gene was found
				#if this line contains the ENST that is for this gene,
				#record this id in found
				for this_id in enst_id:
					if(this_id) in quant_line:
						found=this_id
				#if any(this_id in quant_line for this_id in enst_id):

				#if this id has been found, remove it from the list and
				#write this to output
				if(found !=-1 ):
					print(quant_line)
					#rna_output.write(found)
					enst_id.remove(found) #remove to save time
					rna_output.write(quant_line)
					#if all enst_id are found, stop searching
					if(len(enst_id)==0):
						break				
	#return how many enst ids were not found
	return(len(enst_id))




#########################################################################
outfile=sys.argv[1]
ensg_id=sys.argv[2]
R1=sys.argv[3]
R2=sys.argv[4]
outfolder=sys.argv[5]
enst=sys.argv[6]

#get correct file name for check if file quantification has been performed
if(str(R2)==str("-1")):
	potential_file_name=str("salmon_"+R1)
else:
	potential_file_name=str("salmon_"+R1+"_"+R2)
is_file=str("salmon/"+potential_file_name+"/quant.sf")

#check if salmon quantification has already been performed
#if it has not, quantify
if not os.path.exists(is_file):
	rna_file=do_salmon(R1,R2)
else:
	rna_file=is_file
	print("RNA FILE : "+rna_file)

#get ensts via biomart
enst_id_dic=ensg_to_enst(ensg_id)
enst_ids=list(enst_id_dic)
#if optional enst supplied, add to the list
if not(enst == "-1"):
	enst_ids.append(enst)


#get the enst quantifications and output as a new more informative file (outfolder/outfile)
num_not_found=get_rna_quant(enst_ids,rna_file,outfolder,outfile)

#if none of the enst were found, then print the error
if(num_not_found==len(list(enst_id_dic))):
	print("TRANSCRIPTS NOT FOUND --")
	print(list(enst_id_dic))





