#Many SCMs are identified near repetitive regions which are likely false positive events due to the nature of genomic deletions.
#Reyka Jayasinghe (reykajayasinghe@gmail.com)
#Last Edit: December 8th, 2018

#Identified short introns computationally from this paper:  TableS1 - 
#>HG38_CHR1_1257252_1257277_+ +  >HG38_CHR1_1257245_1257251_+ +  >HG38_CHR1_1257278_1257284_+ +
#CTGGGAGGTTTGAAAGGAAATTCT        ATAGAT  GGGAAA  3
#>HG38_CHR1_1299905_1299973_+ +  >HG38_CHR1_1299898_1299904_+ +  >HG38_CHR1_1299974_1299980_+ +
#CTGTTGGGGGTGGCATTAGGGAAGGTCACGGAGGCCTGGCCCAGCCCCACCCCTCCCAAAGGCTCAC     GGTCCA  CGGCTG  3

#NOTE: Data may need to be lifted over.

#After liftover check below coordinates of SCM with above repetitive sequence coordinates
#chr1:201178736-201178736
#chr1:201178836-201178836
#chr1:201179142-201179142
#chr16:1254349-1254349

#Modules
import sys
import os
from collections import defaultdict
import subprocess
import re
import time

USI="TableS1_USI_repeats.txt"

file=sys.argv[1]
#Determine current working directory

#usidata
usidict=defaultdict(dict)
USI_file=open(USI,"r")
lines = USI_file.readlines()
for i in range(0, len(lines)):
	line = lines[i]
	if line.startswith(">"):
		(usirange,leftflank,rightflank)=line.strip().split('\t')
		usidict[usirange]=lines[i+1].strip()

#liftedover to hg38 mutation data
myfile=open(file,"r")
for line2 in myfile:
	(chrom,crange)=line2.strip().split(':')#chr1:201178736-201178736
	(mutation,stop)=crange.split('-')
	for rusi in usidict:
		(buildusi,chromusi,startusi,endusi,strand)=rusi.split('_')#>HG38_CHR1_1257252_1257277_+	
		mstartusi=int(startusi)-1000
		mendusi=int(endusi)+1000
		if chromusi == chrom.upper():
			if int(mutation) > int(mstartusi) and int(mutation) < int(mendusi):
				print(line2.strip(),rusi,usidict[rusi])
		
#Output example
#PositionOfInterest	OverlappingUSI	Sequence	LeftSeq	RightSeq	Score(Results of Deletion Screener (0 indicates direct repeat, 1 indicates inverted repeat, 2 indicates homopolymeric expansion, and 3 indicates no repetitive sequences were detected))
#chr1:111414946-111414946 >HG38_CHR1_111414895_111414920_+ + CACAGGGGTCACAGACTGATGACC	GATGAC	CAGGGG	0
#chr1:111414946-111414946 >HG38_CHR1_111414902_111414972_+ + GTCACAGACTGATGACCCACAGGGGTCAGGGTCTTTTCCCCAGGGGTCACAGACTGATAACCCACAGAG	CACAGG	CAGGGT	3
#chr1:240207698-240207698 >HG38_CHR1_240207639_240207673_+ + CTCCTCCGCCCCCTCTACCCGGAGCGGCAATAC	GGAATA	CCTCCG	0
#chr22:29489578-29489578 >HG38_CHR22_29489593_29489618_+ + GCCAAGTCCCCAGAGAAGGAAGAG	TGAGAA	CAAGTC	0
#chr22:29489629-29489629 >HG38_CHR22_29489593_29489618_+ + GCCAAGTCCCCAGAGAAGGAAGAG	TGAGAA	CAAGTC	0
			
