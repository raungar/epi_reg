#!/bin/python

import os
import sys
import urllib
import biomartpy as bm


def do_salmon(R1,R2,outfolder):
	#create quants directory for salmon quantification if it odes not exist
	if not os.path.exists(str(outfolder+"/quants")):
		os.mkdir(str(outfolder+"/quants"))

	#perform salmon quantification on single end or paired end data
	print("salmon quanitification")
	if(int(R2)==-1):
		print("single end salmon")
		file_name=str("salmon_"+R1)
		os.system(str("salmon quant -i hs.grch38.index -l A -r "+outfolder+"/rna/"+R1+".fastq.gz -p 8 -o "+outfolder+"/quants/salmon_"+R1))
	else:
		file_name=str("salmon_"+R1+"_"+R2)
		print("paired end salmon")
		os.system(str("salmon quant -i hs.grch38.index -l A -1 "+outfolder+"/rna/"+R1+".fastq.gz -2 "+outfolder+"/rna/"+R2+".fastq.gz -p 8 -o "+outfolder+"/quants/salmon_"+R1+"_"+R2))
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


#NEED A BETTER WAY OF GETTING THIS INFORMATION
def get_rna_quant(enst_id,rna_file,outfolder,outfile):
	#get number of max matches so once all ids are found they can stop running through them
	my_file=str(outfolder+"/quants/"+rna_file+"/quant.sf")

	with open(outfile,"w") as rna_output:
		with open(my_file, "r") as quant_file:
			for quant_line in quant_file:
				found=-1
				for this_id in enst_id:
					if(this_id) in quant_line:
						found=this_id
				#if any(this_id in quant_line for this_id in enst_id):
				if(found !=-1 ):
					print(quant_line)
					#rna_output.write(found)
					enst_id.remove(found) #remove to save time
					rna_output.write(quant_line)
					#if all enst_id are found, stop searching
					if(len(enst_id)==0):
						break				
#def rna_main(ensg_id):
#rna_main("ENSG00000118689")


outfile=sys.argv[1]
ensg_id=sys.argv[2]
R1=sys.argv[3]
R2=sys.argv[4]
outfolder=sys.argv[5]


print("OUTFLIE: "+outfile)

rna_file=do_salmon(R1,R2,outfolder)
enst_id_dic=ensg_to_enst(ensg_id)
get_rna_quant(list(enst_id_dic),rna_file,outfolder,outfile)


