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

#make bins for chip peak analysis
def make_bins(chr,start_init, end_init):
	start=str(start_init-50000)
	end=str(end_init+50000)
	
	os.system(str("echo "+chr)+"\t"+start+"\t"+end+str(" | sed 's/\s/\t/g' > output/chr_file_"+chr+"_"+start+"_"+end+".txt"))
	os.system(str("bedtools makewindows -b output/chr_file_"+chr+"_"+start+"_"+end+'''.txt -w 500 | awk '{fix=$2+1; print $1"\t"fix"\t"$3e}' > output/bin_file_'''+chr+"_"+start+"_"+end+".bed"))
	os.remove("output/chr_file_"+chr+"_"+start+"_"+end+".txt")


#get peak score
def get_peaks(chr,start_init, end_init,chip_file_id,outfile):
	#region is 50,000kb +/- of start and stop
	start=str(start_init-50000)
	end=str(end_init+50000)

	#AGAIN, MAYBE JUST MAKE ALL THESE AT ONCE BEFORE LOOP
	#if the directory does not exist, create it
	if not os.path.exists("output/chip_peaks"):
		os.mkdir("output/chip_peaks")
	
	#get peaks at intersection with bin file
	os.system(str("bedtools intersect -wao -a  output/bin_file_"+chr+"_"+start+"_"+end+".bed -b output/rna/"+chip_file_id+".bed.gz > output/chip_peaks/"+chip_file_id+".bed"))
	#sort for merge step
	os.system(str("sort -k1,1 -k2,2n output/chip_peaks/"+chip_file_id+".bed > output/chip_peaks/sorted_"+chip_file_id+".bed"))
	os.system(str("bedtools merge -i output/chip_peaks/sorted_"+chip_file_id+".bed -c 14 -o max > output/chip_peaks/"+outfile+".bed"))
	os.remove(str("output/chip_peaks/sorted_"+chip_file_id+".bed"))
	os.remove(str("output/chip_peaks/"+chip_file_id+".bed"))



def hg19_to_GRCh38(id):
	if not os.path.exists("output/hg19ToHg38.over.chain.gz"):
		urllib.request.urlretrieve("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", "output/hg19ToHg38.over.chain.gz")
		urllib.request.urlretrieve("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver","output/liftOver")
		os.chmod("output/liftOver",777) #give permissions to operate on
	#os.system(str("output/liftOver output/rna/"+file+".bed.gz output/hg19ToHg38.over.chain.gz output/rna/lifted_"+file+".bed.gz unMapped"))
	print("converting hg19 to grch38")
	read_file=str("output/rna/"+id+".bed.gz")
	write_file=str("output/rna/GRCh38_"+id+".bed.gz")
	lo=LiftOver("output/hg19ToHg38.over.chain.gz")


	with gzip.open(read_file,"rt") as bed_file:
		with gzip.open(write_file,"wt") as write_bed:
			for bedline in bed_file:
				this_bedline=str(bedline).strip().split("\t")
				#grab converted coordinates
				coord1=lo.convert_coordinate(str(this_bedline[0]),int(this_bedline[1]))
				coord2=lo.convert_coordinate(str(this_bedline[0]),int(this_bedline[2]))
				#if not found, just don't print
				if not coord1:
					continue
				if not coord2:
					continue
				#replace coordinates are write to new file
				output_line=("\t".join([str(coord1[0][0]),str(coord1[0][1]),str(coord2[0][1]),this_bedline[3],
						this_bedline[4],this_bedline[5],this_bedline[6],this_bedline[7],
						this_bedline[8],this_bedline[9]]))
				write_bed.write(output_line+"\n")
	print("Done: "+write_file)


#def main():

outfile=sys.argv[1]
bed_id=sys.argv[2]
chr=sys.argv[3]
start=sys.argv[4]
end=sys.argv[5]
assemble=sys.argv[6]

if(assemble=="hg19"):
	hg19_to_GRCh38(bed_id)
	bed_id=str("GRCh38_"+bed_id)

	chr="chr6"
#	start=108559823
#	end=108684774
#	#bed_file="ENCFF998UVS" #ENCFF327LZT
	#bed_file="ENCFF327LZT"
get_peaks(chr, int(start),int(end),bed_id,outfile)

#main()
