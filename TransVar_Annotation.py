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
import time
start = time.time()

file=sys.argv[1]

#Determine current working directory
dir_path = os.path.dirname(os.path.realpath(file))

#Load canonical transcripts for each gene
#File below derived from uniprot_isoform_ids and contains all transcripts for all genes. 
#uniprot_isoform_id	uniprot_entry_id	is_canonical	ensembl_gene_id	ensembl_tx_id	ensembl_protein_id	tx_row_num	tx_rank	hugo_gene_id	symbol	chromosome	start	end	tsl	tx_length_bp
#A0A075B6H9-1	A0A075B6H9	1	ENSG00000211637	ENST00000390282	ENSP00000374817	1	1	5921	IGLV4-69	22	22030934	22031472		419
#A0A075B6I0-1	A0A075B6I0	1	ENSG00000211638	ENST00000390283	ENSP00000374818	1	1	5931	IGLV8-61	22	22098700	22099212		414
#A0A075B6I1-1	A0A075B6I1	1	ENSG00000211639	ENST00000390284	ENSP00000374819	1	1	5920	IGLV4-60	22	22162199	22162681		362

canonical="/Users/rjayasin/Desktop/OneDrive/Programs/GitHub/MiSplice_Supplemental/all_uniprot_ensembl_id_mappings.2018-10-16.tsv"

canontrans=defaultdict(dict)
canonsize=defaultdict(dict)
cfile=open(canonical,"r")
for line in cfile:
	(uniprot_isoform_id,uniprot_entry_id,is_canonical,ensembl_gene_id,ensembl_tx_id,ensembl_protein_id,tx_row_num,tx_rank,hugo_gene_id,symbol,chromosome,start,end,biotype,tsl,tx_length_bp)=line.strip().split('\t')
	#If gene already has a transcript stored in the dictionary, then check if the current transcript size is larger than the one already stored
	if ((tsl == "1") and (biotype == "protein_coding")):
		if symbol in canonsize:
			current_tx_size=canonsize[symbol]
			if int(tx_length_bp) > int(current_tx_size):
				canontrans[symbol]=ensembl_tx_id
				canonsize[symbol]=tx_length_bp
			else:
				continue
		#If no transcript current stored for gene of interest then add it
		else:
			canontrans[symbol]=ensembl_tx_id
			canonsize[symbol]=tx_length_bp
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
					#If canonical transcript is not specified then take first transcript if it is specified as (protein_coding) after transcript name
					else:
						#replace tabs in entry with commas
						alt=alltranscripts[i].split('\t')
						finalt=",".join(alt)
						#append remaining transcript data together
						othertranscriptdata.append(finalt)
		
			#Combine all alt transcripts into one variable
			try:
				finalcanon=canonicaldata[0]
				note = "longest"
			except IndexError:
				note="other"
				firsttranscript=othertranscriptdata[0]
				dataother=firsttranscript.split(',')
				#check to see if transcirpt is (protein_coding)
				currenttrans=dataother[1].split(' ')
				if currenttrans[1] == "(protein_coding)":
					CSQNother=dataother[6].split(';')
					mutclassification=CSQNother[0]
					finalcanon=",".join(dataother) #extract first transcript info for reporting
				else:
					note="Not Reported"
			allalttranscripts="|".join(othertranscriptdata)
			annotated=line.strip()+"\t"+note+"\t"+mutclassification+"\t"+finalcanon+"\t"+allalttranscripts+"\n"
			output.write(annotated)
		
		else:
			print("!!!!!!!ERROR!!!!!!!!!",line)
			errorline="!!!!!!!ERROR!!!!!!!!!\t"+line+"\n"
			errorout.write(errorline)
output.close()
pfile.close()
errorout.close()

end = time.time()
print(end - start)
