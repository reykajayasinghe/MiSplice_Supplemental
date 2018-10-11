#Reyka Jayasinghe (reyka@wustl.edu)
#Last edited: October 11th, 2018
#This script takes in an intermediate output of misplice and does the following:
#1) Filter out specific genes
#2) Combines samples that have the same mutation in the same cancer type into one line entry (https://github.com/ding-lab/misplice)
#3) Filters out SCM events that have > 5% of controls having at least one read with the same SCM event
#4) Requires a minimum of 20 controls
#5) Combines mutations that are linked to the same SCM event. These are put into a separate category: CANCER.rgSCM.multiplemutations
#6) Annotates samples with genomic context information and splice score and saved to: CANCER.rgSCM.filtered
#USAGE: python /gscmnt/gc2706/dinglab/medseq/LabCode/Reyka/MiSplice_Supplemental/Final_Filtering.py ../ACC/novel.splice.scores.rc.key.combined.noHLA.vaf.highexp ACC

#Example file input data
#high_expression  brca	TCGA-TT_T  HLA-A   6_29912108_G_C  304	 0,0,6,0,0,4335,0,0,3847,0,1732,2,0,0,0,0,0,6,0,0,0,1018,2236,0,2,16,0,0,4,1,0,0,1
#TCGA-PK-A5H8_19_53116969_G_A	TCGA-PK-A5H8	19	53116969	G	A	3ss	-	53116956	AATCCACACTGGAGAGAAACCTT	AATCCATACTGGAGAGAAACCTT	-22.36	-22.57	207	TCGA-PK-A5H8_19_53116969_G_A	high_expression	acc	TCGA-PK-A5H8_T	ZNF83	19_53116969_G_A	8	0,0,0,0,0,0,0,0,0,0,0,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
#TCGA-OR-XXXX_8_143958154_G_A	6.65107	TCGA-OR-XXXX_8_143958154_G_A	TCGA-OR-XXXX	8	143958154	G	A	3ss	-	143958144	CTCGCTGGACCAGCCCCAAGGTG	CTCGCTGGATCAGCCCCAAGGTG	-3.93	-2.45	6841	TCGA-OR-XXXX_8_143958154_G_A	high_expression	acc	TCGA-OR-XXXX_T	CYP11B1	8_143958154_G_A	455	0,0,8,0,23,56,0,0,0,0,0,3,0,3,3,0,0,0,0,3,0,0,6,202,60,0,0,2,0,0,19,6,13,0,0,0,3,0,0,109,17,7,0,0,4,0,0,0,16,10,77,0,0,0,0,0,14,0,12,33,0,0,0,33,0,0,14,13,0,0,0,11,0,5,40,24,0,5


#Modules
import sys
import os
from collections import defaultdict

cancertype=sys.argv[2]
file=sys.argv[1]
print("Processing:",file)

#Determine current working directory
dir_path = os.path.dirname(os.path.realpath(file))

print(dir_path)

#Dictionary intitialization
samplecasedict=defaultdict(list)
samplelist=defaultdict(list)
genedict=defaultdict(dict)
controldict=defaultdict(dict)
contextdict=defaultdict(dict)
vafdict=defaultdict(dict)
readdict=defaultdict(dict)

############################
########GENE FILTER#########
############################
#After manually reviewing many potential splice-site-creating mutations in the below genes, many were false positives. Many of the genes int he below list are highly homologous and have very large exons which is a potential reason for the false positive calls.

genes_to_filter=['MUC4','AHNAK','AHNAK2','CRIPAK','IGHV1-17','IGHV1-3','IGHV1-46','IGHV1-69','IGHV3-11','IGHV3-13','IGHV3-23','IGHV3-33','IGHV3-53','IGHV3-64','IGHV3-66','IGHV3-7','IGHV3-73','IGHV3-9','IGHV4-28','IGHV4-4','IGHV4-59','IGHV4-61','MUC17','MUC20','MUC1','MUC13','MUC16','MUC2','MUC5B','MUC21','MUC6','MUC7']

fr=dir_path+"/temp"
output = open(fr,'w')
pfile=open(file,"r")
for line in pfile:
	#(expression,cancer,sample,gene,site,casereads,controlreads)=line.strip().split('\t')
	(key,vaf,key2,sn1,chrom,pos,ref,alt,sstype,something,position,beforemut,aftermut,beforescore,afterscore,totalreads,key3,expression,cancer,sample,gene,site,casereads,controlreads)=line.strip().split('\t')
	if gene not in genes_to_filter:
		samplecasedict[site].append(casereads)
		genedict[site]=gene
		samplelist[site].append(sample)
		controldict[site]=controlreads
		contextinfo=sstype+"\t"+something+"\t"+position+"\t"+beforemut+"\t"+aftermut+"\t"+beforescore+"\t"+afterscore+"\t"
		contextdict[key]=contextinfo
		readdict[key]=totalreads
		vafdict[key]=vaf
	else:
		continue
####################################
#COMBINING MUTATIONS BY CANCER TYPE#
####################################
for sites,values in genedict.items():
	casevalues=','.join(map(str,samplecasedict[sites]))
	casevallist=samplecasedict[sites]
	casesamples=','.join(map(str,samplelist[sites]))
	genename=genedict[sites]
	controlvalues=controldict[sites]
	controllist=controlvalues.split(",")   
	
	readcount=0
	controlcounter=0
#######################################
#######TOTAL CONTROLS FILTER###########
#######################################
#Filter out samples with less than 20 supporting controls
#######################################
#######CONTROL FRACTION FILTER#########
#######################################
#Remove sites with greater than 5% of control samples having the same reported alternative splicing event
#Determine total number of control samples that have SCM event reported without mutation
	for control in controllist:
		controlcounter += 1
		controlint=int(control)
		if controlint > 1:
			readcount += 1
	finalcount=(float(readcount)/float(controlcounter))*100
#Filter 
	if finalcount<5.00:
		if controlcounter > 19:
			printme=str(finalcount)+"\t"+genename+"\t"+sites+"\t"+casesamples+"\t"+cancertype+"\t"+casevalues+"\t"+controlvalues+"\n"
			output.write(printme)

pfile.close()
output.close()
###########################################
##COMBINE NEARBY MUTATIONS INTO ONE ENTRY##
###########################################

tempfile=open(fr,"r")

mutations = defaultdict(dict)
#Example input data
#Fraction_of_Controls_with_Event Gene Mutation Sample Cancer Case_Reads Empty Control_Reads
#0	ZNF93	19_20045120_C_T	TCGA-XX_T	ESCA	29		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
#0	ZNF93	19_20045119_C_A	TCGA-XX_T	ESCA	29		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
#0	ZNF93	19_20045126_G_A	TCGA-YY_T	OV	20		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
#0	ZNF93	19_20044900_T_A	TCGA-YY_T	OV	20		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0

for line in tempfile:
	(controlreadfraction,gene,position,sample,cancer,reads,controlreads)=line.strip().split("\t")
	currentsamplelist=sample.split(",")
	currentsamplelist.sort()
	key=gene+"_"+str(currentsamplelist)
	#If gene and sample list already exist then there is a nearby mutation for this cancer-sample set that we need to combine
	if key in mutations:
		currentdata=mutations[key]
		(controlreadfraction1,gene1,position1,sample1,cancer1,reads1,controlreads1)=currentdata.strip().split("\t")
		#Check 1 - make sure sample list is the same
		samplelist=sample1.split(",")
		samplelist.sort()
		readslist=reads1.split(",")
		if set(currentsamplelist) == set(samplelist):
			pass
		else:
			continue
		#Check 2 - make sure total reads are the same
		if set(reads) == set(reads1):
			pass
		else:
			continue
		#Check 3 - combine and replace
		newline=controlreadfraction+","+controlreadfraction1+"\t"+gene+"\t"+position+","+position1+"\t"+sample+"\t"+cancer+"\t"+reads+"\t"+controlreads
		mutations[key]=newline
	#Gene sample pair does not already exist
	else:
		mutations[key]=line.strip()
tempfile.close()
#######################
#####OUTPUT DATA#######
#######################

out = cancer+".rgSCM.filtered"
out2 = cancer+".rgSCM.multiplemutations.filtered"
output = open(out,'w')
output2 = open(out2,'w')

header="Key\tFraction_of_Controls_with_Event\tGene\tMutation\tTotal_Mutations\tSample_List\tTotal_Samples\tCancer\tCase_Reads\tTotal_Case_Depth\tJAF\tType_Splice_Site\tStrand\tSCM_Genomic_Position\tBefore_Mutation_Context\tAfter_Mutation_Context\tBefore_Score\tAfter_Score\tControl_Reads\n"
#OUTPUT
#Key	Fraction_of_Controls_with_Event	Gene	Mutation	Total_Mutations	Sample_List	Total_Samples	Cancer	Case_Reads	Control_Reads
#EIF2AK1_['TCGA-4Z-AA7R_T']	0.0	EIF2AK1	7_6094337_T_C	1	TCGA-4Z-AA7R	1	BLCA	83	231	35.9307	3ss	-	6094318	AGAATCTGATGTTCCAGCAGAAA	GGAATCTGATGTTCCAGCAGAAA	3.86	3.53		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
output.write(header)

for key,value in mutations.items():
	(crf,g,pos,sn,can,read,controlre)=value.strip().split("\t")
	#Count number of samples for this mutation
	finalsamplelist=sn.split(",")
	slistlength=len(finalsamplelist)
	#Count number of mutations for this SCM
	poslist=pos.split(",")
	poslistlength=len(poslist)
	#Empty out variables for each iteration
	datacontext=""
	datavaf=""
	final=""
	vafdata=""
	keysingle=""
	#Extract other context info (if only one sample and one mutation, context should be same)
	if poslistlength == 1:
		if slistlength == 1:
			sn = sn[:-2] #Remove TCGA-XX-XXXX_T (remove _T) to compare sample names
			keysingle=sn+"_"+pos
			datacontextval=contextdict[keysingle]
			datavaf=vafdict[keysingle]
			datareads=readdict[keysingle]
			#Final print statement
			final=key+"\t"+crf+"\t"+g+"\t"+pos+"\t"+str(poslistlength)+"\t"+sn+"\t"+str(slistlength)+"\t"+can+"\t"+read+"\t"+str(datareads)+"\t"+str(datavaf)+"\t"+str(datacontextval)+"\t"+controlre+"\n"   
			output.write(final)
	#If more than one sample but the same mutation, then only the vaf should differ.
		if slistlength > 1:
			vafs=[]
			keynew=""
			reads=[]
			for s in finalsamplelist:
				#redefine new key for each sample to extract appropriate VAF info
				s = s[:-2]
				keynew=s+"_"+pos
				vafdata=vafdict[keynew]
				vafs.append(vafdict[keynew])
				reads.append(readdict[keynew])
			datacontext=contextdict[keynew]
			readsvar=','.join(reads)
			vafsvar=','.join(vafs)
			final=key+"\t"+crf+"\t"+g+"\t"+pos+"\t"+str(poslistlength)+"\t"+sn+"\t"+str(slistlength)+"\t"+can+"\t"+read+"\t"+readsvar+"\t"+vafsvar+"\t"+str(datacontext)+"\t"+controlre+"\n"		
			output.write(final)
	#If there is more than one mutation, context needs to be recomputed.
	if poslistlength > 1:
		final=key+"\t"+crf+"\t"+g+"\t"+pos+"\t"+str(poslistlength)+"\t"+sn+"\t"+str(slistlength)+"\t"+can+"\t"+read+"\t"+controlre+"\n"   			
		output2.write(final)
							

output.close()
output2.close()

print("Completed, check final output CANCER.rgSCM.filtered")

######################################################
########ADD BACK IN MUTATION ANNOTATION INFO##########
######################################################
#Use TransVar to annotate all variants.http://bioinformatics.mdanderson.org/main/Transvar

##MAF info was lost in these samples, add back in for this final list.


#Cleanup
os.remove(fr)
