#Annotate MiSplice post-filtered results with TransVar, Adds two columns to the end of the input file format that includes the canonical transcript results and all alternative transcript results.
#Reyka Jayasinghe (reykajayasinghe@gmail.com)
#Last Edit: October 15th, 2018

##Download TransVar
#https://transvar.readthedocs.io/en/latest/download_and_install.html
#sudo pip install transvar #download transvar
#transvar config --download_ref --refversion hg19 #Download reference
#transvar config --download_anno --refversion hg19 #Set up databases
#transvar config -k reference -v [path_to_hg19.fa] --refversion hg19 #link reference to transvar if you already have one

###Annotation of Genomic Coordinates
#https://transvar.readthedocs.io/en/latest/annotation_from_genomic_level.html

#input file format:
#Key	Fraction_of_Controls_with_Event	Gene	Mutation	Total_Mutations	Sample_List	Total_Samples	Cancer	Case_Reads	Total_Case_Depth	JAF	Type_Splice_Site	Strand	SCM_Genomic_Position	Before_Mutation_Context	After_Mutation_Context	Before_Score	After_Score	Control_Reads
#A2ML1_['TCGA-XX_T']	0.982800982801	A2ML1	12_8975820_C_T	1	TCGA-XX	1	BLCA	5	35	14.2857	3ss	+	8975826	TAAATTTCCCCTCCGTTCAGAAG	TAAATTTCCCCTCTGTTCAGAAG	6.66	6.70	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
#AAAS_['TCGA-XX_T']	1.0989010989	AAAS	12_53715207_G_T	1	TCGA-XX	1	ESCA	7	58	12.069	5ss	-	53715209	GGGGTCAAG	GGGGTAAAG	-6.89	3.53	0,0,0,2,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0

#Modules
import sys
import os
from collections import defaultdict
import subprocess
import re

file=sys.argv[1]

#Determine current working directory
dir_path = os.path.dirname(os.path.realpath(file))

#Load canonical transcripts for each gene
#ZYX:ENST00000322764
#ZZEF1:ENST00000381638
#ZZZ3:ENST00000370801
canonical="/Users/rjayasin/Desktop/OneDrive/Programs/GitHub/MiSplice_Supplemental/E75_bed_v3_canonical_transcripts.tsv"
canontrans=defaultdict(dict)
cfile=open(canonical,"r")
for line in cfile:
	(gene,trans)=line.strip().split(':')
	canontrans[gene]=trans
cfile.close()

ofile=file+".transvar"
output=open(ofile,"w")
errorout=open("error.transvar","w")

pfile=open(file,"r")

for line in pfile:
	if not line.startswith("Key"):
		(key,fraction,gene,mutation,total_mutations,sample_list,total_samples,cancer,case_reads,total_case_depth,JAF,Type_Splice_Site,Strand,SCM_gposition,beforemutation,aftermutation,beforescore,afterscore,empty,control_reads)=line.strip().split('\t')
		#For transvar genomic coordinate annotation
		#SNV : transvar ganno --ensembl -i 'chr3:g.178936091G>A'
		#Insertion: transvar ganno -i 'chr2:g.69741762_69741763insTGC' --ensembl
		#Inframe Deletion : transvar ganno -i "chr2:g.234183368_234183379del" --ensembl
		#Frameshift Deletion : transvar ganno -i "chr2:g.234183368_234183380del" --ensembl
		#Putting genomic cooridnates in format for transvar command line 
		mutation=mutation.split('_')
		refsize=len(mutation[2])
		mutsize=len(mutation[3])
		valid_DNA="ACGT"
		gDNA=""
		#SNP
		if (refsize == 1 and mutsize == 1):
			if (mutation[2] in valid_DNA and mutation[3] in valid_DNA):
				gDNA="chr"+mutation[0]+":g."+mutation[1]+mutation[2]+">"+mutation[3]
				print("SNP:",gDNA)
		#Insertion
		if (mutation[2] == "-" and mutsize >= 1):
			start=int(mutation[1])
			end=start+1
			gDNA="chr"+mutation[0]+":g."+str(start)+"_"+str(end)+"ins"+mutation[3]
			print("Insertion:",gDNA)
		#Deletion
		if (refsize >= 1 and mutation[3] == "-"):
			start=int(mutation[1])
			end=start+(refsize-1)
			gDNA="chr"+mutation[0]+":g."+str(start)+"_"+str(end)+"del"
			print("Deletion:",gDNA)
		if (gDNA is not ''):
			result = subprocess.run(['transvar', 'ganno', '-i', gDNA, '--ensembl'], stdout=subprocess.PIPE)
			#print(result.stdout)
			alltranscripts=result.stdout.decode('utf-8').split('\n')
			#total number of transcripts
			size=len(alltranscripts)
			i=0
			#Remove last empty line element from transcript list
			alltranscripts.pop()
			#Empty variables/lists for storing transcript data
			canonicaldata=[]
			othertranscriptdata=[]
			mutclassification=""
			#Extract the canonical transcript info and concatenate all other transcript info into one entry
			for i in range(len(alltranscripts)):
				if i == 0: #skip header
					continue
				else:
					(input,transcript,gene,strand,coordinates,region,info)=alltranscripts[i].strip().split('\t')
					trans=transcript.strip().split(' ')
					canontranscript=canontrans[gene]
					if canontranscript==trans[0]:
						#Extract CSQN (consequence)
						data=alltranscripts[i].split('\t')
						CSQN=data[6].split(';')
						mutclassification=CSQN[0]	
						fincanon=",".join(data)
						canonicaldata.append(fincanon)
					else:
						#replace tabs in entry with commas
						alt=alltranscripts[i].split('\t')
						finalt=",".join(alt)
						#append remaining transcript data together
						othertranscriptdata.append(finalt)
		
			#Combine all alt transcripts into one variable
			try:
				finalcanon=canonicaldata[0]
			except IndexError:
				finalcanon="Not Reported"
			allalttranscripts="|".join(othertranscriptdata)
			annotated=line.strip()+"\t"+mutclassification+"\t"+finalcanon+"\t"+allalttranscripts+"\n"
			output.write(annotated)
		
		else:
			print("!!!!!!!ERROR!!!!!!!!!",line)
			errorline="!!!!!!!ERROR!!!!!!!!!\t"+line+"\n"
			errorout.write(errorline)
output.close()
pfile.close()
errorout.close()
