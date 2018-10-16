
# MiSplice Supplemental Scripts

Reyka Jayasinghe (reyka@wustl.edu)
Last Edit: October 15th, 2018

### Final_Filtering.py
This script takes in an intermediate output of misplice(https://github.com/ding-lab/misplice) and does the following:
1) Filter out subset of genes
2) Combines samples that have the same mutation in the same cancer type into one line entry
3) Filters out SCM events that have > 5% of controls having at least one read with the same SCM event
4) Requires a minimum of 20 controls
5) Combines mutations that are linked to the same SCM event. These are put into a separate file: CANCER.rgSCM.multiplemutations
6) Annotates samples with genomic context information and splice score and saved to: CANCER.rgSCM.filtered.txt
```
USAGE: python Final_Filtering.py novel.splice.scores.rc.key.combined.noHLA.vaf.highexp ACC
```

### TransVar Annotation
Annotate MiSplice post-filtered results with TransVar. Adds two columns to the end of the input file format that includes the canonical transcript results and all alternative transcript results.

#### Download TransVar
```
https://transvar.readthedocs.io/en/latest/download_and_install.html
sudo pip install transvar #download transvar
transvar config --download_ref --refversion hg19 #Download reference
transvar config --download_anno --refversion hg19 #Set up databases
transvar config -k reference -v [path_to_hg19.fa] --refversion hg19 #link reference to transvar if you already have one
```
#### Annotation of Genomic Coordinates
-Error output file: error.transvar
-Final output: CANCER.rgSCM.filtered.txt.transvar

```
USAGE: python TransVar_Annotation.py CANCER.rgSCM.filtered
```

