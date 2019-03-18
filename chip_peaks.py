#!/bin/python

import os
import sys
import stat
import gzip

#this file creates the chip peak quantifications in directory outfolder/chip_peaks

#get peak score
def get_peaks(chr,start_init, end_init,chip_file_id,outfile,outfolder):
	#region is 50,000kb +/- of start and stop
	start=str(start_init-50000)
	end=str(end_init+50000)

	#if the directory does not exist, create it
	if not os.path.exists(str(outfolder+"/chip_peaks")):
		os.mkdir(str(outfolder+"/chip_peaks"))

	#get peaks at intersection with bin file
	os.system(str("bedtools intersect -wao -a  "+outfolder+"/bin_file_"+chr+"_"+start+"_"+end+".bed -b chip/"+chip_file_id+".bed.gz > "+outfolder+"/chip_peaks/intersect_"+chip_file_id+".bed"))
	#sort for merge step
	os.system(str("sort -k1,1 -k2,2n "+outfolder+"/chip_peaks/intersect_"+chip_file_id+".bed > "+outfolder+"/chip_peaks/sorted_"+chip_file_id+".bed"))
	#perform merge and output into final file
	os.system(str("bedtools merge -i "+outfolder+"/chip_peaks/sorted_"+chip_file_id+".bed -c 14 -o max > "+outfile+".bed"))
	print(str("MERGED: "+outfile+".bed"))
	#remove the temporary files
	os.remove(str(outfolder+"/chip_peaks/sorted_"+chip_file_id+".bed"))
	os.remove(str(outfolder+"/chip_peaks/intersect_"+chip_file_id+".bed"))


#get input
outfile=sys.argv[1]
bed_id=sys.argv[2]
chr=sys.argv[3]
start=sys.argv[4]
end=sys.argv[5]
assemble=sys.argv[6]
outfolder=sys.argv[7]


#skip if hg19
if(assemble=="hg19"):
	print("NOT FUNCTIONAL CURRENTLY")
else:
	#get the peaks and print
	get_peaks(chr, int(start),int(end),bed_id,outfile,outfolder)
