#Reyka Jayasinghe (reyka@wustl.edu)
#Last edited: October 15th, 2018

#Key	PercentControlswithASEvent	Gene	Position	Total Mutations	Sample List	Total Samples with Mutation	Cancer Type	Case Reads	Total Case Depth	Junction Allele Fraction	Type Splice Site	Strand	Site	Distance from Novel to Canonical Site	Canonical Site Position	Canonical Site Pre Mutation	Canonical Site Post Mutation	Canonical Site Score Pre Mutation	Canonical Site Score Post Mutation	Novel Site Position	Novel Site Pre Mutation	Novel Site Post Mutation	Novel Site Score Pre Mutation	Novel Site Score Post Mutation	Population Frequency	Controls	AnnotationNotes	Consequence	Annotated Transcript	All Transcripts 
#C11orf83_['TCGA-TT-TTTT_T']	4.31893687708	C11orf83	11_62439424_G_A	1	TCGA-TT-TTTT	1	CESC	17	238	7.14286	3ss	+	TCGA-TT-TTTT_11_62439424_G_A	3	62439424	ACTCCCCTTCCTGCCCTCAGGAG	ACTCCCCTTCCTGCCCTCAAGAG	10.74	1.99	62439427	CCCCTTCCTGCCCTCAGGAGATG	CCCCTTCCTGCCCTCAAGAGATG	-7.45	-4.12	0.000352609	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,67,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,2,0,0,0,3,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0	other	CSQN=IntronicSNV	chr11:g.62439424G>A,ENST00000528862 (protein_coding),C11ORF48,-,chr11:g.62439424G>A/c.1-4283C>T/.,inside_[5-UTR;intron_between_exon_1_and_2],CSQN=IntronicSNV;aliases=ENSP00000434489;source=Ensembl	chr11:g.62439424G>A,ENST00000528862 (protein_coding),C11ORF48,-,chr11:g.62439424G>A/c.1-4283C>T/.,inside_[5-UTR;intron_between_exon_1_and_2],CSQN=IntronicSNV;aliases=ENSP00000434489;source=Ensembl|chr11:g.62439424G>A,ENST00000531323 (protein_coding),C11ORF83,+,chr11:g.62439424G>A/c.121-1G>A/.,inside_[intron_between_exon_2_and_3],CSQN=SpliceAcceptorSNV;C2=SpliceAcceptorOfExon3_At_chr11:62439424;aliases=ENSP00000432692;source=Ensembl|chr11:g.62439424G>A,ENST00000377953 (protein_coding),C11ORF83,+,chr11:g.62439424G>A/c.121-1G>A/.,inside_[intron_between_exon_1_and_2],CSQN=SpliceAcceptorSNV;C2=SpliceAcceptorOfExon2_At_chr11:62439424;aliases=ENSP00000367189;source=Ensembl

#Modules
import sys
import os
import subprocess
from collections import defaultdict
from subprocess import call
from subprocess import Popen, PIPE
from subprocess import PIPE, run


#Dictionaries
Mutationinfo=defaultdict(dict)
beddata=defaultdict(dict)
bedinfo=defaultdict(dict)
genomicpositions=defaultdict(list)
CDSsize=defaultdict(int)
exoncount=defaultdict(int)
Chromlist=[]

#Hard coded files
#hg19 reference sequence
hg19genome="/Users/rjayasin/Desktop/hs37d5.fa"

#Functions

def sizecalc( str ):
	#passed string: chr1_46740389_46743487
	try:
		(c,s,e,strand)=str.split('_')
		size=abs(int(e)-int(s)+1)
		return size
	except AttributeError: #this is usually the case when the mutation is in the last exon so there is no downstream intron to grab
		return "no data"

def mutation_distance ( mutationposition,exoninfo ):
	(c,start,stop,strand)=exoninfo.split('_')
	distostart=int(mutationposition)-int(start)
	distoend=int(mutationposition)-int(end)
	if strand == "+":
		return distostart,distoend
	if strand == "-":
		return distoend,distostart

def eidata ( eiinfo ):
	#Define upstream and downstream intron
	#-------------*------------ mut position
	#e1----i1----e2----i2----e3 gene structure (+) 
	#e3----i2----e2----i1----e1 gene structure (-)
	#(+ and - strand) if mut is in e2, grab i1(upstreamintron) and i2(downstreamintron)
	exonnum=eiinfo.replace("e","")
	idown=exonnum
	iup=int(exonnum)-1
	upstream_intron=gene+":"+transcript+":"+strand+":i"+str(iup)
	downstream_intron=gene+":"+transcript+":"+strand+":i"+str(idown)
	upintron=bedinfo.get(upstream_intron)
	downintron=bedinfo.get(downstream_intron)
	return upintron,downintron,exonnum

def CDS_length ( gene ):
	positions=genomicpositions.get(gene)
	length=0
	for exon in positions:
		chrom,start,stop,strand=exon.split("_")
		current_exonsize=abs(int(stop)-int(start))
		length+=current_exonsize
	return length

def reverse_complement(dna):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '|': '|'}
	#return ''.join([complement[base] for base in dna]) #complement
	return ''.join([complement[base] for base in dna[::-1]]) #reverse_complement

def grab_additional_sequence ( startnew,endnew,strand ):
#if region extends into intron, we need to grab additional bp sequence
	mutgenomicrange=chrom.replace("chr","")+":"+str(startnew)+"-"+str(endnew)
	proc2 = subprocess.Popen(['samtools', 'faidx', hg19genome, mutgenomicrange], stdout=subprocess.PIPE)
	output2 = proc2.stdout.read().decode("utf-8")
	nucleotides_intron=output2.split("\n")
	nucleotides_intron.pop(0)
	mutexonnew=[]
	if strand == "-":
		nucleotide_intron_reversed=nucleotides_intron[::-1]
		for seqi in nucleotide_intron_reversed:
			mutexonnew.insert(0,seqi)
		mutseqnew="".join(mutexonnew)
	if strand ==  "+":
		for seqi in nucleotides_intron:
			mutexonnew.append(seqi)
		mutseqnew="".join(mutexonnew)
	return mutseqnew	


def grab_sequence (gene,region_bedformat,mutposition,typess,newsspos,TOI,Refallele,Altallele):
	positions=genomicpositions.get(gene)
	WTsequence=[]
	MUTsequence=[]

	for exon in positions:
		(chrom,start,stop,strand)=exon.split("_")
		#Defined genomic range is different based on if derived from positive or negative strand
		#Need to modify start of gene coding sequence based on MET start site
		#Currently bed file contains 5' and 3' UTR info
		#Determine where start codon is located
		if TOI in CDSstart:
			startsite,cdsstartsite=CDSstart[TOI].split(",") #Key:ENST00000371975 Value:46714077,46713403
		else: #No start site reported in database! Skip this an an annotation error.
			startsite="unknown"
			break
		#Modify start of first exon to have true transcriptional start
		if strand == "+" and int(startsite) < int(stop) and int(startsite) > int(start):#if true start site is before the end of the current exon, modify current exon start position
			start=startsite
		if strand == "+" and int(startsite) > int(stop) and int(startsite) > int(start):#if true start site is located downstream of current exon, skip this exon cause it is only UTR
			continue
		if strand == "-" and int(startsite) < int(stop) and int(startsite) > int(start):#if true start site is before the end of the current exon, modify current exon start position
			stop=int(startsite)-1
		if strand == "-" and int(startsite) < int(stop) and int(startsite) < int(start):#if true start site is located downstream of current exon, skip this exon cause it is only UTR
			continue
		
		else:
			pass
		
		modifiedstart=int(start)+1
		modifiedstop=int(stop)+1
		genomicrange=chrom.replace("chr","")+":"+str(modifiedstart)+"-"+str(modifiedstop)
		#Add or remove "chr" depending on input fasta sequence
		#Example output
		#>1:46743489-46743652
		#GCAGGGACCATTGAGGAGAAGATCTTCCAGCGTCAGAGCCACAAGAAGGCACTGAGCAGC
		#TGTGTGGTGGATGAGGAGCAGGATGTAGAGCGCCACTTCTCTCTGGGCGAGTTGAAGGAG
		#CTGTTTATCCTGGATGAAGCTAGCCTCAGTGACACACATGACAG
		#call(command, shell=True)
		proc = subprocess.Popen(['samtools', 'faidx', hg19genome, genomicrange], stdout=subprocess.PIPE)
		output = proc.stdout.read().decode("utf-8")
		nucleotides=output.split("\n") #['>1:46739810-46739888', 'GTACTTATACGTCCGCCTGGATGGCACGATGTCCATTAAGAAGCGAGCCAAGGTTGTAGA', 'ACGCTTCAATAGTCCATCG', '']
		#Need to reverse order in which we add DNA bases together because sequence is negative strand derived
		if strand == "-":	
			nucleotide_order=nucleotides[::-1]
		#determine exon with mutation, this is where the sequence will need to be altered for everything downstream
		if (nucleotides[0] == region_bedformat):
			mutexon=[]
			if strand == "+":
				for seq in nucleotides:
					mutexon.append(seq)
					#WTsequence.append(seq)
				bpposition=mutexon.pop(0).split(":")
				bpposition.pop(0)
				(exonstart,exonend)=bpposition[0].split("-")
				lengthseq=abs(int(exonend)-int(exonstart))
				if typess == "5ss":
					distfromcanon=abs(int(newsspos)-int(exonstart)+1)
					distfromend=int(exonend)-int(newsspos)
					wtseq="".join(mutexon)
					if int(distfromend) > 0: #exon variant  
						mutseq=wtseq[0:distfromcanon]
					if int(distfromend) < 0: #intron variant
						#if region extends into intron, we need to grab additional bp sequence
						newend=int(newsspos)
						mutseq=grab_additional_sequence(exonstart,newend,strand)
					#If mutation is contained within new exon we need to modify mutsequence to contain mutation
					if int(mutposition) < int(newsspos):
						mutseq_alt=list(mutseq) #split string into list and replace Ref with Alt allele
						distfrommuttostart=int(mutposition)-int(exonstart)
						refcheck=mutseq_alt[distfrommuttostart]
						if refcheck == Refallele:
							mutseq_alt[distfrommuttostart]=Altallele
							mutseq_m="".join(mutseq_alt)
						else:
							print(distfrommuttostart,distfromcanon,mutposition,newsspos,refcheck,Refallele,"error wtf!")
							sys.exit()
						mutseq=mutseq_m
					if int(mutposition) > int(newsspos):
						pass
					newexonsize=distfromcanon
					MUTsequence.append(mutseq)
					MUTsequence.append("|")
					WTsequence.append(wtseq)
					WTsequence.append("|")	
				if typess == "3ss":
					distfromcanon=int(newsspos)-int(exonstart)+2
					distfrommuttostart=int(mutposition)-int(exonstart)
					wtseq="".join(mutexon)
					newsspos_actualstartexon=int(newsspos)+1
					if int(distfromcanon) > 0: #exon variant
						mutseq=wtseq[distfromcanon:]
					if int(distfromcanon) < 0: #intron variant
						newstart=int(newsspos)+2
						mutseq=grab_additional_sequence(newstart,modifiedstop,strand)
					if int(mutposition) > int(newsspos_actualstartexon): #within new exon
						distfromnewss=int(mutposition)-int(newsspos)-2
						mutseq_alt=list(mutseq)
						refcheck=mutseq_alt[distfromnewss]
						if refcheck == Refallele:
							mutseq_alt[distfromnewss]=Altallele
							mutseq_m="".join(mutseq_alt)
						else:
							print(distfrommut,distfromcanon,mutposition,newsspos,refcheck,Refallele,"error wtf!")
							sys.exit()
						mutseq=mutseq_m
					if int(mutposition) < int(newsspos):
						pass
					distfromend=abs(int(newsspos)-int(exonend)+1)
					newexonsize=distfromend
					MUTsequence.append(mutseq)
					MUTsequence.append("|")
					WTsequence.append(wtseq)
					WTsequence.append("|")	
			if strand == "-":	
				for seq in nucleotide_order:
					#For Mutant sequence we need to determine where the new splice site occurs and shift the coding sequence
					mutexon.insert(0,seq)
					WTsequence.insert(0,seq)
				bpposition=mutexon.pop(0).split(":")
				WTsequence.pop(0) #remove >1:46743489-46743652 from sequence
				bpposition.pop(0) #remove chr : >1
				(exonstart,exonend)=bpposition[0].split("-")
				lengthseq=abs(int(exonend)-int(exonstart))
				#For negative strand start will be a smaller number than start
				# ==============---------------=============
				#               3ss(AG)       5ss(GT)  
				if typess == "5ss": #chr19	37681052	37697940
					distfromcanon=int(exonstart)-int(newsspos)-2
					wtseq="".join(mutexon)
					if int(distfromcanon) < 0: # exon variant
						mutseq=wtseq[abs(distfromcanon):]
					if int(distfromcanon) > 0: # intron variant
						#if region extends into intron, we need to grab additional bp sequence
						newstart=int(newsspos)+2
						mutseq=grab_additional_sequence(newstart,modifiedstop,strand)
					#Mutation is contained within new exon need to modify mutseq to contain mutation
					if int(mutposition) > int(newsspos):
						distfrommut=int(mutposition)-int(newsspos)-2
						mutseq_alt=list(mutseq) #split string into list and replace Ref with Alt allele
						refcheck=mutseq_alt[distfrommut]
						if refcheck == Refallele:
							mutseq_alt[distfrommut]=Altallele
							mutseq_m="".join(mutseq_alt)
						else:
							print(distfrommut,distfromcanon,mutposition,newsspos,refcheck,Refallele,"error wtf!")
							sys.exit()
						mutseq=mutseq_m
					if int(mutposition) < int(newsspos):
						pass	
					newexonsize=abs(int(exonend)-int(newsspos)-1)
					#mutseq=wtseq[distfromcanon:]
					MUTsequence.insert(0,mutseq)
					MUTsequence.insert(0,"|")
				if typess == "3ss":
					distfromend=int(exonend)-int(newsspos)-1
					distfromcanon=int(exonstart)-int(newsspos)-1
					wtseq="".join(mutexon)
					if int(distfromend) < 0: #intron variant
						#newend=int(newsspos)+2	
						#mutseq=grab_additional_sequence(modifiedstart,newend,strand)	
						mutseq=grab_additional_sequence(modifiedstart,newsspos,strand)
					if int(distfromend) > 0: #exon variant
						#mutseq=wtseq[0:distfromcanon]
						newend2=int(abs(distfromcanon))
						#mutseq=wtseq[0:newsspos]
						mutseq=wtseq[0:newend2]
					#If Mutation is contained within new exon need to modify mutseq to contain mutation
					if int(mutposition) < int(newsspos):
						distfrommut=int(mutposition)-int(exonstart)
						mutseq_alt=list(mutseq)
						refcheck=mutseq_alt[distfrommut]
						if refcheck == Refallele:
							mutseq_alt[distfrommut]=Altallele
							mutseq_m="".join(mutseq_alt)
						else:
							print(mutposition,newsspos,"error wtf!")
							sys.exit()
						#print("B:",mutseq)
						#print("A:",mutseq_m)
						mutseq=mutseq_m
					if int(mutposition) > int(newsspos):
						pass
					#print("MUT:",typess,distfromcanon,distfromend,exonstart,newsspos,newend2,mutseq)
					#print("WT:",typess,distfromcanon,distfromend,exonstart,newsspos,newend2,wtseq)
					newexonsize=abs(int(exonstart)-int(newsspos)-1)
					#MUTsequence.append(mutseq) #append reversed sequence so it is in positive strand orientation
					MUTsequence.insert(0,mutseq)
					MUTsequence.insert(0,"|")
		else:
			if strand == "-":
				#just concatenate sequence and add to final sequence		
				nucleotide_order.pop(0)#Remove first element '>1:46739810-46739888'
				nucleotide_order.pop(-1)#Remove last element (empty) , '']
				#nucleotides.insert(0,"|") #this will mark the end of the exon designating the exon-exon junction for NMD analysis later on
				otherexon=[]
				for seq in nucleotide_order:
					otherexon.insert(0,seq)
					#WTsequence.insert(0,seq)
					#MUTsequence.insert(0,seq)
				exonseq="".join(otherexon)
				MUTsequence.insert(0,exonseq)
				MUTsequence.insert(0,"|")
				WTsequence.insert(0,exonseq)
				WTsequence.insert(0,"|")
			if strand == "+":
				nucleotides.pop(0)#Remove first element '>1:46739810-46739888'
				nucleotides.pop(-1)#Remove last element (empty) , '']
				nucleotides.append("|") #this will mark the end of the exon designating the exon-exon junction for NMD analysis later on
				#otherexon=[]
				for seq in nucleotides:
					WTsequence.append(seq)
					MUTsequence.append(seq)
	#Return mutant and wildtype sequences
	ifinalmutseq="".join(MUTsequence)
	ifinalwtseq="".join(WTsequence)
	if strand == "-":
		finalwtseq=reverse_complement(ifinalwtseq)	
		finalmutseq=reverse_complement(ifinalmutseq)
	if strand == "+":
		finalwtseq=ifinalwtseq
		finalmutseq=ifinalmutseq
	try:	
		finalmutseq = finalmutseq[:-1] #remove last | from end of sequence
		finalwtseq = finalwtseq[:-1]
	except UnboundLocalError:
		finalmutseq = "Annotation Error"
		finalwtseq = "Annotation Error"
		newexonsize = "Annotation Error"
	try:
		#print("WT:",finalwtseq)
		#print("MUTT:",finalmutseq)
		return finalmutseq,finalwtseq,newexonsize
	except UnboundLocalError:
		if (startsite =="unknown"):
			return "No Reported Start","No Reported Start","No Reported Start"
		else:
			return "5UTRvariant","5UTRvariant","5UTRvariant" #Sometimes transvar annotates mutaitons of the donor or acceptor site that are upstream as exon 1 as Splice donor or acceptor variants when they should be classified as UTR variants. We do not analyze UTR variants in this study so these will be skipped for further analyses.
#Amino acid dictionary table, adapted from: https://pythonforbiologists.com/dictionaries/
gencode = { 'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 
 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
 'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
 'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',	
 'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
 'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
 'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}



def translate ( sequence,TOI ):
	totalbp=len(sequence)
	AAsequence=[]
	#Determine where start codon is located
	#if TOI in CDSstart:
	#	startsite,cdsstartsite=CDSstart[TOI].split(",") #Key:ENST00000371975 Value:46714077,46713403
	#else:
#		startsite="unknown"
	#starting from start codon, translate sequence
#	firststartcodon=abs(int(cdsstartsite)-int(startsite))
	#determine location of all exon junctions
	EJ="|"
	ejpos=[pos for pos, char in enumerate(sequence) if char == EJ] #[677, 765, 886, 948, 1085, 1156, 1446, 1572, 1724, 1852, 1928, 2060, 2172, 2297, 2377, 2558, 2665, 3067]
	seqwithoutej=sequence.replace("|","")#remove | from sequence for translation	

	#split sequence by every third character
	#lookup every three characters in gencode table for amino acid position	
	for i in range(0,totalbp,3):
		j=i+3
		bp=seqwithoutej[i:j]
		if bp in gencode:
			AA=gencode[bp]
			AAsequence.append(AA)
			if AA == "*": #once a stop codon is reached break out of translation
				ptcpos=i
				break
	translatedsequence="".join(AAsequence)
	######################################
	####NONSENSE MEDIATED DECAY CHECK#####
	######################################
	#If a premature termination codon resides >= 50 nt upstream of an exon-exon junction
	#then NMD will successfully degrade the aberrant transcript. 
	#determine where stop codon is introduced
	#if PTC is located > 50 nucleotides upstream of last exon-exon junction --> should elicit NMD
	#If PTC is in last exon or within 50 nucleotides of last exon-exon junction --> escape NMD degradation
	if len(ejpos) > 0:

		lastej=ejpos[-1]
	
		try:
			if int(ptcpos) > int(lastej):
				disttoej=abs(int(lastej)-int(ptcpos))
				NMD="ESCAPE-NMD"
			if int(ptcpos) < int(lastej):
				#check if within 50 bp of last exon-exon junction
				disttoej=abs(int(lastej)-int(ptcpos))
				if int(disttoej) < 50:
					NMD="ESCAPE-NMD"
				if int(disttoej) >= 50:
					NMD="ELICIT-NMD"
			if int(ptcpos) == int(lastej):
				NMD="ESCAPE-NMD"
		#If mutation induces readthrough of stop codon due to frameshift event
		except UnboundLocalError:
			NMD="ESCAPE-NMD"
		return translatedsequence,NMD		
	#This will be the case for genes with only one exon Eg. HIST1H2BK
	if len(ejpos) == 0:
		NMD="ESCAPE-NMD"
		return translatedsequence,NMD		

#####################################
#STORE EXON and INTRON INFO FROM CDS#
#####################################
bedfile="E75_bed_v3.longest.aav2.tsv" #longest transcript
enst_transcripts=[]#List of ENSEMBL transcripts of interest
#INPUT BED FILE
#chr1    2005157 2005367 PRKCZ:ENST00000495347:+:e1
#chr1    2005368 2029784 PRKCZ:ENST00000495347:+:i1
#chr1    2029785 2029993 PRKCZ:ENST00000495347:+:e2
#chr1    2029994 2066699 PRKCZ:ENST00000495347:+:i2
bed=open(bedfile,"r")
for line in bed:
	(chrom,start,stop,geneinfo)=line.strip().split('\t')
	(genename,transcript,strand,info)=geneinfo.split(":")
	enst_transcripts.append(transcript)
	key=chrom+"_"+start+"_"+stop+"_"+strand
	beddata[key]=geneinfo	
	bedinfo[geneinfo]=key
	if info.startswith("e"):
		exon_num=int(info.replace("e",""))
		exonsize=sizecalc(key)
		#Tracking total number of exons for gene of interest
		if genename in exoncount:
			val=exoncount.get(genename)
			if int(exon_num) > int(val):
				exoncount[genename]=exon_num
		else:
			exoncount[genename]=exon_num
		#Adding all exon coordinates to list
		genomicpositions[genename].append(key)	
bed.close()
#get unique transcripts
unique_enst=set(enst_transcripts)

###########################################
#STORE UN TRANSLATED REGION INFO FROM UCSC#
###########################################
#Table Browser; https://genome.ucsc.edu/cgi-bin/hgTables
#clade: Mammal, Genome: Human, assembly: Feb. 2009 (GRCh37/hg19)
#group: gene and gene predictions, track: GENCODE gene V27lift37
#table: basic, outputformat: BED, getoutput > Whole gene
utrbed="RefSeqUCSC_genes_GRCh37.tsv"

#chrom	txStart	txEnd	transcript	bin	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
#chr1	46713403	46744144	ENST00000371975.8_1	0	+	46714077	46743954	0	18	677,87,120,61,136,70,289,125,151,127,75,131,111,124,79,180,164,401,	0,780,2268,10954,12232,12810,12995,13529,19727,22927,24734,24940,25623,25892,26406,26806,30085,30340,

CDSstart=defaultdict(dict)

utr=open(utrbed,"r")
for line in utr:
	data=line.strip().split('\t')
	tscript=data[3].split('_')#ENST00000371975.8_1	 
	tscriptonly=tscript[0].split('.')
	if tscriptonly[0] in unique_enst:
		if data[5] == "+":
			CDSstart[tscriptonly[0]]=data[6]+","+data[1] #Key:ENST00000371975 Value:46714077,46713403
		if data[5] == "-":
			CDSstart[tscriptonly[0]]=data[7]+","+data[1] #Key:ENST00000371975 Value:46714077,46713403
utr.close()

#################################################
####EXTRACT REGION OF INTEREST FROM MUT FILE#####
#################################################
header="Key\tMutation\tEffect\tEffectSize\tFrame\tUpstreamIntronSize\tDownstreamIntronSize\tExonSize\tNewExonSize\tExonFraction\tDistancetoStart\tDistancetoEnd\tCDSLength\tTotalExons\tNMDPrediction\tMutantAA\tWTAA\tKey\tPercentControlswithASEvent\tGene\tPosition\tTotal Mutations\tSample List\tTotal Samples with Mutation\tCancer Type\tCase Reads\tTotal Case Depth\tJunction Allele Fraction\tType Splice Site\tStrand\tSite\tDistance from Novel to Canonical Site\tCanonical Site Position\tCanonical Site Pre Mutation\tCanonical Site Post Mutation\tCanonical Site Score Pre Mutation\tCanonical Site Score Post Mutation\tNovel Site Position\tNovel Site Pre Mutation\tNovel Site Post Mutation\tNovel Site Score Pre Mutation\tNovel Site Score Post Mutation\tPopulation Frequency\tControls\tAnnotationNotes Consequence\tAnnotated Transcript\tAll Transcripts\n"
file=sys.argv[1]
dat=open(file,"r")
err=file+".error.sites"
out=file+".translated.sites"
ERRFILE=open(err,"w")
OUTPUT=open(out,"w")
OUTPUT.write(header)

num_lines = sum(1 for line in open(file))

current_count=0

for line in dat:
	current_count+=1
	data=line.strip().split('\t')
	newsplicesite=int(data[20])-1 #is this true for both positive and negative strand sites?
	originalsplicesite=int(data[15])-1
	typesplicesite=data[11]
	(chrom,mutposition,ref,alt)=data[3].split('_') #data[3] = 1_46743596_C_T
	chrom_chr="chr"+chrom+"_"
	Mutationinfo[data[13]]=line.strip() #data[4] = TCGA-TN-APPP_1_42345596_C_T
	(chromosomed,p,ref,alt)=data[3].split('_') #1_42345596_C_T
	if chrom not in Chromlist:
		Chromlist.append(chrom)
	#go through bed dictionary and determine where mutation is located
	for region in beddata:
		if region.startswith(chrom_chr):
			(ch,start,end,strand)=region.split("_")
			#Determine if mutation falls within exon or intron region
			#########EXON SHRINKAGE ONLY#############
			expandedstart=int(start)-10
			expandedend=int(end)+10
			#if (int(mutposition) >= int(start) and int(mutposition) <= int(end)):
			if (int(mutposition) >= int(expandedstart) and int(mutposition) <= int(expandedend)):
				value=beddata.get(region)
				#Determine if mutation is in intron or exon
				(gene,transcript,strand,eiinfo)=value.split(':') #RAD54L:ENST00000371975:+:e17
				#print(value,region,eiinfo,expandedstart,mutposition,expandedend)
				#if eiinfo.startswith('i'):
				#	continue
				if not eiinfo.startswith('i'):
					(upintron,downintron,currentexonnum)=eidata(eiinfo)
					upsize=sizecalc(upintron)
					downsize=sizecalc(downintron)
					esize=sizecalc(region)
					cdslength=str(CDS_length(gene))
					totalexons=exoncount.get(gene)	
					#exon fraction = current exon/totalexons
					exonfraction=(int(currentexonnum)/int(totalexons))*100
					(disttostart,disttoend)=mutation_distance(mutposition,region)
					#Convert region (chr1_46743488_46743651_+) to compare with bed file (>1:46743489-46743652)
					region_bedformat=ch.replace("chr",">")+":"+str(int(start)+1)+"-"+str(int(end)+1)
					(mutseq,wtseq,newexonsize)=grab_sequence(gene,region_bedformat,mutposition,typesplicesite,newsplicesite,transcript,ref,alt)
					if mutseq=="5UTRvariant":
						err="Annotation Issue Skipping: 5UTR variant"+"\t"+data[0]+"\t"+data[3]+"\n"
						ERRFILE.write(err)
						break
					if mutseq=="Annotation Error":
						err="Annotation Issue Skipping: Check Gene Transcript Overlap"+"\t"+data[0]+"\t"+data[3]+"\n"
						ERRFILE.write(err)	
						break
					if mutseq=="No Reported Start":
						err="Annotation Issue Skipping: No Reported CDS Start Site"+"\t"+data[0]+"\t"+data[3]+"\n"
						ERRFILE.write(err)
						break
					mutaatranslation,mutnmd=translate(mutseq,transcript)	
					wtaatranslation,wtnmd=translate(wtseq,transcript)
					difference=int(esize)-int(newexonsize)
					effect=""
					effectsize=abs(difference)
					if int(difference) > 0:
						effect="exon-shrinkage"
					if int(difference) < 0:
						effect="exon-extension"
					frame=""
					calc=difference/3

					if int(difference) % 3 == 0:
						frame="in-frame"
					else:
						frame="off-frame"
					printstatement=data[0]+"\t"+data[3]+"\t"+effect+"\t"+str(effectsize)+"\t"+frame+"\t"+str(upsize)+"\t"+str(downsize)+"\t"+str(esize)+"\t"+str(newexonsize)+"\t"+str(exonfraction)+"\t"+str(disttostart)+"\t"+str(disttoend)+"\t"+str(cdslength)+"\t"+str(totalexons)+"\t"+mutnmd+"\t"+mutaatranslation+"\t"+wtaatranslation+"\t"+line.strip()+"\n"
					OUTPUT.write(printstatement)
	
	print(data[0],str(current_count),"/",num_lines)
dat.close()
ERRFILE.close()
OUTPUT.close()

#### Working draft of scirpt to load in your own fasta file and store whole genome. Under development.
#genome=open(hg19genome,'r')
#hg19=defaultdict(dict)
#for line in genome:
#	line.strip()
#	#genomelist=list(line)
#	if line.startswith('>'): #>1 dna:chromosome chromosome:GRCh37:1:1:249250621:1
#		info=line.strip().split(' ')
#		chrom=info[0] #>1
#		chrnum=chrom[1:] #1
#	else:
#		if chrnum in uniquechrlist:
#			if chrnum in hg19:
#				currentval=hg19[chrnum]
#				newval=currentval+line.upper()
#				hg19[chrnum]=newval
#			else:
#				hg19[chrnum]=line.upper()	
#		else:
#			continue
#print("Finished storing hg19 genome")



