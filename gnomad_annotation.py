#Reyka Jayasinghe (reyka@wustl.edu)
#Edited: October 16th, 2018
#Usage: python gnomad_annotation.py mutationfile
#Note: Only works for snps, if indels are in mutation file it will throw an error

#To download gnomad data
#http://gnomad.broadinstitute.org/downloads
#mkdir gnomad_data

#Grab following linnes from gnomad VCF files for annotation to Mutation File
##INFO=<ID=AC_POPMAX,Number=A,Type=Integer,Description="AC in the population with the max AF">
##INFO=<ID=AF_AFR,Number=A,Type=Float,Description="Allele Frequency among African/African American genotypes, for each ALT allele, in the same order as li
##INFO=<ID=AF_AMR,Number=A,Type=Float,Description="Allele Frequency among Admixed American genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF_ASJ,Number=A,Type=Float,Description="Allele Frequency among Ashkenazi Jewish genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF_EAS,Number=A,Type=Float,Description="Allele Frequency among East Asian genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF_FIN,Number=A,Type=Float,Description="Allele Frequency among Finnish genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF_NFE,Number=A,Type=Float,Description="Allele Frequency among Non-Finnish European genotypes, for each ALT allele, in the same order as listed
##INFO=<ID=AF_OTH,Number=A,Type=Float,Description="Allele Frequency among Other (population not assigned) genotypes, for each ALT allele, in the same orde

#Example INPUT file
#mutationfile (make sure this has unique mutations)
#1	100382185	G	A
#1	100741158	G	A
#1	104166597	G	A


#Example OUTPUT
#mutationfile.gnomad
#1.909222.G.A    AF_POPMAX=6.33173e-02   AF_AFR=6.33173e-02      AF_AMR=2.38663e-03      AF_ASJ=0.00000e+00      AF_EAS=0.00000e+00      AF_FIN=0.00000e+00      AF_NFE=6.69344e-05      AF_OTH=1.02041e-03
#1.909253.C.T    AF_POPMAX=2.67523e-04   AF_AFR=0.00000e+00      AF_AMR=0.00000e+00      AF_ASJ=0.00000e+00      AF_EAS=0.00000e+00      AF_FIN=0.00000e+00      AF_NFE=2.67523e-04      AF_OTH=0.00000e+00
#1.985359.C.T    AF_POPMAX=3.81503e-03   AF_AFR=3.81503e-03      AF_AMR=1.19332e-03      AF_ASJ=0.00000e+00      AF_EAS=0.00000e+00      AF_FIN=0.00000e+00      AF_NFE=0.00000e+00      AF_OTH=0.00000e+00


#Import necessary packages
from collections import defaultdict
import sys
import gzip
import time

start = time.time()

#Directory where Gnomad Data is located
#File name scheme: gnomad.genomes.r2.0.2.sites.chr1.vcf.bgz
GNOMAD_DIR="/gscmnt/gc2532/dinglab/gnomad_data/"
GNOMAD_FILENAME="gnomad.genomes.r2.0.2.sites.chr"

#Dictionary
mutations = defaultdict(list)

#Mutation File 
#Format: 1   242048615   G   C
file=sys.argv[1]
MUT=open(file,"r")
f=file+".gnomad"
OUT=open(f,"w")
#Store all chromosome,position,ref,alt allele information in dictionary by chromosome
#Example Key value pair:
#(Key,:,Values (list format))
#('11', ':', ['11.8251891.G.A', '11.111959590.G.T'])
#('10', ':', ['10.72358655.G.A', '10.72358655.G.A', '10.72358655.G.A', '10.72358655.G.A'])

for line in MUT:
	(chromosome,position,refallele,altallele) = line.strip().split('\t')
#Skip all indels
	if len(refallele) > 1:
		continue
	if len(altallele) > 1:
		continue
	key = chromosome+"."+position+"."+refallele+"."+altallele
	mutations[chromosome].append(key)
MUT.close()
for chrom,sites in mutations.items():
	#Get unique set of sites
	sites_unique=set(sites)
		
	#open associated gnomad VCF file to find allele frequency information
	VCF=GNOMAD_DIR+GNOMAD_FILENAME+chrom+".vcf.bgz"
	print(chrom,":",sites,",",VCF)
	#FILE=open(VCF)
	with gzip.open(VCF,'r') as fin:
		for vline in fin:
			if not vline.startswith("#"):
				data=vline.strip('\'').split('\t')
				vcf_key=data[0]+"."+data[1]+"."+data[3]+"."+data[4]
				if vcf_key in sites_unique:
					data2=vline.split(';')
					#Since each VCF line is different we need to go through each element in line and store all the pop allele frequencies based on name.
					for entry in data2:
						if entry.startswith("AF_AFR"):
							AF_AFR=entry
						if entry.startswith("AF_AMR"):
							AF_AMR=entry
						if entry.startswith("AF_ASJ"):
							AF_ASJ=entry
						if entry.startswith("AF_EAS"):
							AF_EAS=entry
						if entry.startswith("AF_FIN"):
							AF_FIN=entry
						if entry.startswith("AF_NFE"):
							AF_NFE=entry
						if entry.startswith("AF_OTH"):
							AF_OTH=entry
						if entry.startswith("AF_POPMAX"):
							AF_POPMAX=entry
							popmax=AF_POPMAX.split('=')
#							frequency=(float(popmax[1])*100)	
						if entry.startswith("AF_FilterStatus"):
							AS_FilterStatus=entry
						
					thing=vcf_key+"\t"+AF_POPMAX+"\t"+AF_AFR+"\t"+AF_AMR+"\t"+AF_ASJ+"\t"+AF_EAS+"\t"+AF_FIN+"\t"+AF_NFE+"\t"+AF_OTH+"\n"
					OUT.write(thing)
				else:
					continue
	fin.close()

OUT.close()

end = time.time()
print(end - start)
