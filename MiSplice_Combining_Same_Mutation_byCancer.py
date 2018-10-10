#Reyka Jayasinghe (reyka@wustl.edu)
#Last edited: October 10th, 2018
#This script takes in an intermediate output of misplice and combines samples that have the same mutation in the same cancer type into one line entry (https://github.com/ding-lab/misplice)

#Example file input data
#high_expression  brca    TCGA-TT_T  HLA-A   6_29912108_G_C  304     0,0,6,0,0,4335,0,0,3847,0,1732,2,0,0,0,0,0,6,0,0,0,1018,2236,0,2,16,0,0,4,1,0,0,1

#Modules
import sys
from collections import defaultdict

file=sys.argv[1]
print("Processing:",file)
#file="novel.splice.scores.rc.key.combined.vaf.noHLA.highexp.maf.casecontrolonly" #intermediate output from https://github.com/ding-lab/misplice
#file="ACC.novel.splice.scores.rc.key.combined.vaf.noHLA.highexp.maf.casecontrolonly"

#Dictionary intitialization
samplecasedict=defaultdict(list)
samplelist=defaultdict(list)
genedict=defaultdict(dict)
controldict=defaultdict(dict)

fr=file+".addfraction"
output = open(fr,'w')

pfile=open(file,"r")
for line in pfile:
    (expression,cancer,sample,gene,site,casereads,controlreads)=line.strip().split('\t')
    samplecasedict[site].append(casereads)
    genedict[site]=gene
    samplelist[site].append(sample)
    controldict[site]=controlreads

for sites,values in genedict.items():
    casevalues=','.join(map(str,samplecasedict[sites]))
    casevallist=samplecasedict[sites]
    casesamples=','.join(map(str,samplelist[sites]))
    genename=genedict[sites]
    controlvalues=controldict[sites]
    controllist=controlvalues.split(",")   
    
    readcount=0
    controlcounter=0

#Determine total number of control samples that have SCM event reported without mutation
    for control in controllist:
        controlcounter += 1
        controlint=int(control)
        if controlint > 1:
            readcount += 1
    finalcount=(float(readcount)/float(controlcounter))*100
    printme=str(finalcount)+"\t"+genename+"\t"+sites+"\t"+casesamples+"\t"+casevalues+"\t"+controlvalues+"\n"
    output.write(printme)

pfile.close()
