#!/bin/python

from Bio import SeqIO
import os
import stat
import gzip
import urllib.request
#import filter_metadata
import tarfile
import zipfile

def do_kallisto(s,l):

	current_dir=os.getcwd()
	os.chdir(str(current_dir+"/output"))
	print(str("\n"+"Building index in: "+os.getcwd()))
	#os.system("kallisto index -i transcripts.idx Homo_sapiens.GRCh38.cdna.all.fa.gz")
	print(str("\n"+"Beginning kallisto quantification in: "+os.getcwd()+"\n"))
	os.system(str("kallisto quant -i hs.grch39_robs.index.gz -o quant_output -b 100 --single rna/ENCFF163DLM.fastq.gz -s "+str(s)+"-l "+l))
	#os.system(str("kallisto quant -i transcripts.idx -o quant_output -b 100 --single rna/ENCFF163DLM.fastq.gz -s "+str(s)+"-l "+l))
	print("kallisto finished")
	os.chdir(str(current_dir))

def do_salmon(paired,R1,R2):

	if not os.path.exists("output/hs.grch39.index"):
		print("building salmon index")
		os.system("salmon index -t output/Homo_sapiens.GRCh38.cdna.all.fa.gz -i output/hs.grch39.index")
	if not os.path.exists("output/quants"):
		os.mkdir("output/quants")
	print("salmon quanitification")
	if(paired=="False"):
		print("single end salmon")
		os.system(str("salmon quant -i output/hs.grch39.index -l A -r output/rna/"+R1+".fastq.gz -p 8 -o output/quants/"+R1))
	else:
		os.system(str("salmon quant -i output/hs.grch39.index -l A -1 output/rna/"+R1+".fastq.gz -2 output/rna/"+R2+".fastq.gz -p 8 -o output/quants/"+R1+"_"+R2))


def get_single_end_params():
	#get file
	if not os.path.exists("output/prinseq.tar.gz"):
		url="https://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz"
		urllib.request.urlretrieve(url,"output/prinseq.tar.gz")
		#	archive = zipfile.ZipFile("output/prinseq.tar.gz", 'r')
		#	imgfile = archive.open('img_01.png')
		t = tarfile.open('output/prinseq.tar.gz', 'r')
		for member in t.getmembers():
			if "prinseq-lite.pl" in member.name:
				t.extract(member, "output")

#make bins for chip peak analysis
def make_bins(chr,start_init, end_init):
	start=str(start_init-50000)
	end=str(end_init+50000)
	
	os.system(str("echo "+chr)+"\t"+start+"\t"+end+str(" | sed 's/\s/\t/g' > output/chr_file_"+chr+"_"+start+"_"+end+".txt"))
	os.system(str("bedtools makewindows -b output/chr_file_"+chr+"_"+start+"_"+end+'''.txt -w 500 | awk '{fix=$2+1; print $1"\t"fix"\t"$3e}' > output/bin_file_'''+chr+"_"+start+"_"+end+".bed"))
	os.remove("output/chr_file_"+chr+"_"+start+"_"+end+".txt")


#get peak score
def get_peaks(chr,start_init, end_init,chip_file_id):
	start=str(start_init-50000)
	end=str(end_init+50000)
	
	#get peaks at intersection with bin file
	os.system(str("bedtools intersect -wao -a  output/bin_file_"+chr+"_"+start+"_"+end+".bed -b output/rna/"+chip_file_id+".bed.gz > output/chip_peaks/"+chip_file_id+".bed"))
	#sort for merge step
	os.system(str("sort -k1,1 -k2,2n output/chip_peaks/"+chip_file_id+".bed > output/chip_peaks/sorted_"+chip_file_id+".bed"))
	os.system(str("bedtools merge -i output/chip_peaks/sorted_"+chip_file_id+".bed -c 14 -o max > output/chip_peaks/merged_"+chip_file_id+".bed"))
	os.remove(str("output/chip_peaks/sorted_"+chip_file_id+".bed"))
	os.remove(str("output/chip_peaks/"+chip_file_id+".bed"))

def main():
	#file_metadata.main()
#	print("HI")
#	with gzip.open("output/rna/ENCFF000JYB.fastq.gz","rt") as file:
#		for record in SeqIO.parse(file, "fastq"):
#		        print(record.id)
#

#	get_single_end_params()
	if not os.path.exists("output/Homo_sapiens.GRCh38.cdna.all.fa.gz"): 
		print("Downloading transcripto fasta")
		#os.system("wget -P output ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz")
		os.system("wget -P output ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz")

		print("Download complete")
	if not os.path.exists("output/chip_peaks"):
		os.mkdir("output/chip_peaks")

	if not os.path.exists("output/hg19ToHg38.over.chain.gz"):
		urllib.request.urlretrieve("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", "output/hg19ToHg38.over.chain.gz")
		urllib.request.urlretrieve("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver","output/liftOver")
		os.chmod("output/liftOver",777) #give permissions

	#s=Estimated average fragment length
	#l=Estimated average fragment length
	#these are kallisto parameters, required for single end.
	#if paired end kallisto calculates this, pass -1, -1 to symbolize paired end
	s=200
	l=10
	#do_kallisto(s,l)
	paired="False"
	#do_salmon(paired, "ENCFF163DLM", -1)
	do_salmon(paired, "ENCFF001RMC,-1)




	chr="chr6"
	start=108559823
	end=108684774
	#bed_file="ENCFF998UVS" #ENCFF327LZT
	bed_file="ENCFF327LZT"
#	make_bins(chr,int(start), int(end)) #do before loop, should only need ot be done once
#	get_peaks(chr, int(start),int(end),bed_file)

main()
