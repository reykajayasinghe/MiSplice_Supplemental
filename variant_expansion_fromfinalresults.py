#For samples sharing the same variant, expand each sample to it's own row. This will be helpful for plotting using R for example.
#Currently only one variant -cancer type pair is represented in the final file, but if you want only one sample-variant-cancer type per row
#you need to use this script to expand the data
#Reyka Jayasinghe (reykajayasinghe@gmail.com)
#Last Edit: December 18ith, 2018

#Modules
import sys
import os
from collections import defaultdict
import subprocess
import re
import time


file=sys.argv[1]

myfile=open(file,"r")
for line in myfile:
	d=line.strip().split('\t')
	samples=d[26].split(',')
	totalsampleswithmutation=len(samples)
	if len(samples)>1:
		creads=d[29].split(',')
		dept=d[30].split(',')
		jafs=d[31].split(',')
		for i in range(len(samples)):
			d[26]=samples[i]
			d[29]=creads[i]
			d[30]=dept[i]
			d[31]=jafs[i]
			print('\t'.join(d))
	else:
		print('\t'.join(d))
