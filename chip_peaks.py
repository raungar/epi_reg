#!/bin/python

import os
import sys
import stat
import urllib
import biomartpy as bm
import subprocess
import pandas as pd
from pyliftover import LiftOver
import gzip


#get peak score
def get_peaks(chr,start_init, end_init,chip_file_id,outfile,outfolder):
	#region is 50,000kb +/- of start and stop
	start=str(start_init-50000)
	end=str(end_init+50000)

	print(str("in get peaks outfolder: "+outfolder+" and outfile: "+outfile))

	#AGAIN, MAYBE JUST MAKE ALL THESE AT ONCE BEFORE LOOP
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





def convert_hg19(id):
	if not os.path.exists("output/hg19ToHg38.over.chain.gz"):
		urllib.request.urlretrieve("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", "output/hg19ToHg38.over.chain.gz")
		urllib.request.urlretrieve("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver","output/liftOver")
		os.chmod("output/liftOver",777) #give permissions to operate on
	#os.system(str("output/liftOver chip/"+file+".bed.gz output/hg19ToHg38.over.chain.gz output/chip/lifted_"+file+".bed.gz unMapped"))
	print("converting hg19 to grch38")
	read_file=str("output/chip/"+id+".bed.gz")
	write_file=str("output/chip/GRCh38_"+id+".bed.gz")
	lo=LiftOver("output/hg19ToHg38.over.chain.gz")
	

	#convert_command=str("output/liftOver <( awk \'{for ( i=1; i<=NF; i++){ if(i==7 || i==8 || i==9){$i=$i*100000} printf $i; if(i==NF){print \"\"} else{printf \"\\t\"}}}\' <( gzip -dc output/chip/"+id+".bed.gz ) ) output/hg19ToHg38.over.chain.gz output/chip/GRCh38_"+id+".bed.gz output/UnMapped")
	#convert_command=str("awk \'{for (i=1; i<=NF; i++){ if(i==7 || i==8 || i==9){$i=$i*100000} printf $i; if(i==NF){print \"\"} else{printf \"\\t\"}}}\' <( gzip -dc output/chip/"+id+".bed.gz ) | head")
	#convert_command = '''bash -c \" awk \'{for (i=1; i<=NF; i++){ if(i==7 || i==8 || i==9){$i=$i*100000} printf $i; if(i==NF){print \"\"} else{printf \"\\t\"}}}\' <( gzip -dc output/chip/ENCFF327LZT.bed.gz ) | head \" '''
	#convert_command = '''awk \'{for (i=1; i<=NF; i++){ if(i==7 || i==8 || i==9){$i=$i*100000} printf $i; if(i==NF){print \"\"} else{printf \"\\t\"}}}\' <( gzip -dc output/chip/ENCFF327LZT.bed.gz ) | head '''
	#convert_command = " bash -c \\" awk \'{for (i=1; i<=NF; i++){ if(i==7 || i==8 || i==9){$i=$i*100000} printf $i; if(i==NF){print \"\"} else{printf \"\\t\"}}}\' <( gzip -dc output/chip/ENCFF327LZT.bed.gz ) | head \\" "
	#awk '{for (i=1; i<=NF; i++){ if(i==7 || i==8 || i==9){$i=$i*100000} printf $i; if(i==NF){print ""} else{printf "\t"}}}' <( gzip -dc output/chip/ENCFF327LZT.bed.gz ) | head
	#convert_command= "awk '{for (i=1; i<=NF; i++){ if(i==7 || i==8 || i==9){\$i=\$i*100000} printf \$i; if(i==NF){print \"\"} else{printf \"\t\"}}}' <( gzip -dc output/chip/ENCFF327LZT.bed.gz ) | head  "
	#convert_command="awk '{for (i=1; i<=NF; i++){ if(i==7 || i==8 || i==9){\$i=\$i*100000} printf \$i; if(i==NF){print \"\"} else{printf \"\t\"}}}' <( gzip -dc output/chip/ENCFF327LZT.bed.gz ) | head  "
	#convert_command = " awk \'{for (i=1; i<=NF; i++){ if(i==7 || i==8 || i==9){\$i=\$i*100000} printf \$i; if(i==NF){print \\"\\"} else{printf \"\\t\"}}}\' <( gzip -dc output/chip/ENCFF327LZT.bed.gz ) | head "



	print("converting")
	#os.system(convert_command)
	#print(str("str=\""+convert_command+"\"; bash -c \"$str\""))
	#subprocess.call(str("str=\""+convert_command+"\"; bash -c \"$str\""), shell=True)
	#print(str("str=\""+convert_command+"\"; bash -c \"$str\""))
	#subprocess.call(str("str=\""+convert_command+"\"; echo \"$str\""), shell=True)






def hg19_to_GRCh38(id):
	if not os.path.exists("output/hg19ToHg38.over.chain.gz"):
		urllib.request.urlretrieve("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", "output/hg19ToHg38.over.chain.gz")
		urllib.request.urlretrieve("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver","output/liftOver")
		os.chmod("output/liftOver",777) #give permissions to operate on
	#os.system(str("output/liftOver output/chip/"+file+".bed.gz output/hg19ToHg38.over.chain.gz output/chip/lifted_"+file+".bed.gz unMapped"))
	print("converting hg19 to grch38")
	read_file=str("output/chip/"+id+".bed.gz")
	write_file=str("output/chip/GRCh38_"+id+".bed.gz")
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
	write_file.close()

#def main():

outfile=sys.argv[1]
bed_id=sys.argv[2]
chr=sys.argv[3]
start=sys.argv[4]
end=sys.argv[5]
assemble=sys.argv[6]
outfolder=sys.argv[7]




print(str("outfile="+outfile+ ", bed_id="+outfile+" , chr="+chr))
print(str("start="+start+" , end="+end+" , assemble=" +assemble +" , outfolder="+outfolder))

if(assemble=="hg19"):
	print("NOT FUNCTIONAL CURRENTLY")
	#convert_hg19(bed_id)
	#bed_id=str("GRCh38_"+bed_id)

else:
	get_peaks(chr, int(start),int(end),bed_id,outfile,outfolder)
	
	chr="chr6"
#	start=108559823
#	end=108684774
#	#bed_file="ENCFF998UVS" #ENCFF327LZT
	#bed_file="ENCFF327LZT"

#main()
