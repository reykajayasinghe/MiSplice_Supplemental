#Reyka Jayasinghe (reyka@wustl.edu)
#Edited: October 16th, 2018
#Usage: python gnomad_annotation.py mutationfile
#Note: Only works for snps, if indels are in mutation file it will throw an error

#To download gnomad data
#http://gnomad.broadinstitute.org/downloads
#mkdir gnomad_data

#Grab following linnes from gnomad VCF files for annotation to Mutation File #2.0 data
##INFO=<ID=AC_POPMAX,Number=A,Type=Integer,Description="AC in the population with the max AF">
##INFO=<ID=AF_AFR,Number=A,Type=Float,Description="Allele Frequency among African/African American genotypes, for each ALT allele, in the same order as li
##INFO=<ID=AF_AMR,Number=A,Type=Float,Description="Allele Frequency among Admixed American genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF_ASJ,Number=A,Type=Float,Description="Allele Frequency among Ashkenazi Jewish genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF_EAS,Number=A,Type=Float,Description="Allele Frequency among East Asian genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF_FIN,Number=A,Type=Float,Description="Allele Frequency among Finnish genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF_NFE,Number=A,Type=Float,Description="Allele Frequency among Non-Finnish European genotypes, for each ALT allele, in the same order as listed
##INFO=<ID=AF_OTH,Number=A,Type=Float,Description="Allele Frequency among Other (population not assigned) genotypes, for each ALT allele, in the same orde

#2.1 data - grabbing all Alternate allele frequency values for each specified population. #Only exome data downloaded, much faster
##INFO=<ID=controls_AF_afr_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of African-American ancestry in the controls subset">
##INFO=<ID=non_neuro_AF_eas_kor,Number=A,Type=Float,Description="Alternate allele frequency in samples of Korean ancestry in the non_neuro subset">
##INFO=<ID=non_topmed_AF_amr,Number=A,Type=Float,Description="Alternate allele frequency in samples of Latino ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_asj_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of Ashkenazi Jewish ancestry in the non_cancer subset">
##INFO=<ID=non_cancer_AF_oth_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of uncertain ancestry in the non_cancer subset">
##INFO=<ID=non_neuro_AF_asj_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of Ashkenazi Jewish ancestry in the non_neuro subset">
##INFO=<ID=non_neuro_AF_afr_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of African-American ancestry in the non_neuro subset">
##INFO=<ID=controls_AF_nfe_swe,Number=A,Type=Float,Description="Alternate allele frequency in samples of Swedish ancestry in the controls subset">
##INFO=<ID=non_neuro_AF_afr_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of African-American ancestry in the non_neuro subset">
##INFO=<ID=non_topmed_AF_amr_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of Latino ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples in the non_cancer subset">
##INFO=<ID=non_cancer_AF_nfe_onf,Number=A,Type=Float,Description="Alternate allele frequency in samples of non-Finnish but otherwise indeterminate European ancestry in the non_cancer subset">
##INFO=<ID=non_cancer_AF_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples in the non_cancer subset">
##INFO=<ID=non_topmed_AF_oth_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of uncertain ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_sas_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of South Asian ancestry in the non_cancer subset">
##INFO=<ID=non_neuro_AF_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples in the non_neuro subset">
##INFO=<ID=controls_AF_afr,Number=A,Type=Float,Description="Alternate allele frequency in samples of African-American ancestry in the controls subset">
##INFO=<ID=non_neuro_AF_eas_jpn,Number=A,Type=Float,Description="Alternate allele frequency in samples of Japanese ancestry in the non_neuro subset">
##INFO=<ID=non_cancer_AF_amr_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of Latino ancestry in the non_cancer subset">
##INFO=<ID=controls_AF_fin_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of Finnish ancestry in the controls subset">
##INFO=<ID=non_neuro_AF_nfe_nwe,Number=A,Type=Float,Description="Alternate allele frequency in samples of North-Western European ancestry in the non_neuro subset">
##INFO=<ID=non_topmed_AF_nfe_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of non-Finnish European ancestry in the non_topmed subset">
##INFO=<ID=non_neuro_AF_sas,Number=A,Type=Float,Description="Alternate allele frequency in samples of South Asian ancestry in the non_neuro subset">
##INFO=<ID=non_cancer_AF_fin_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of Finnish ancestry in the non_cancer subset">
##INFO=<ID=non_cancer_AF_nfe_seu,Number=A,Type=Float,Description="Alternate allele frequency in samples of Southern European ancestry in the non_cancer subset">
##INFO=<ID=non_neuro_AF_nfe_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of non-Finnish European ancestry in the non_neuro subset">
##INFO=<ID=non_neuro_AF_afr,Number=A,Type=Float,Description="Alternate allele frequency in samples of African-American ancestry in the non_neuro subset">
##INFO=<ID=controls_AF_raw,Number=A,Type=Float,Description="Alternate allele frequency in samples in the controls subset, before removing low-confidence genotypes">
##INFO=<ID=non_cancer_AF_eas,Number=A,Type=Float,Description="Alternate allele frequency in samples of East Asian ancestry in the non_cancer subset">
##INFO=<ID=non_cancer_AF_amr_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of Latino ancestry in the non_cancer subset">
##INFO=<ID=non_neuro_AF_nfe_swe,Number=A,Type=Float,Description="Alternate allele frequency in samples of Swedish ancestry in the non_neuro subset">
##INFO=<ID=controls_AF_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples in the controls subset">
##INFO=<ID=non_topmed_AF_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples in the non_topmed subset">
##INFO=<ID=controls_AF_eas_jpn,Number=A,Type=Float,Description="Alternate allele frequency in samples of Japanese ancestry in the controls subset">
##INFO=<ID=controls_AF_nfe_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of non-Finnish European ancestry in the controls subset">
##INFO=<ID=non_neuro_AF_amr,Number=A,Type=Float,Description="Alternate allele frequency in samples of Latino ancestry in the non_neuro subset">
##INFO=<ID=non_neuro_AF_eas_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of East Asian ancestry in the non_neuro subset">
##INFO=<ID=controls_AF_nfe_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of non-Finnish European ancestry in the controls subset">
##INFO=<ID=non_neuro_AF_fin,Number=A,Type=Float,Description="Alternate allele frequency in samples of Finnish ancestry in the non_neuro subset">
##INFO=<ID=non_topmed_AF_sas,Number=A,Type=Float,Description="Alternate allele frequency in samples of South Asian ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_nfe_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of non-Finnish European ancestry in the non_cancer subset">
##INFO=<ID=non_cancer_AF_asj,Number=A,Type=Float,Description="Alternate allele frequency in samples of Ashkenazi Jewish ancestry in the non_cancer subset">
##INFO=<ID=controls_AF_nfe,Number=A,Type=Float,Description="Alternate allele frequency in samples of non-Finnish European ancestry in the controls subset">
##INFO=<ID=controls_AF_oth_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of uncertain ancestry in the controls subset">
##INFO=<ID=controls_AF_asj,Number=A,Type=Float,Description="Alternate allele frequency in samples of Ashkenazi Jewish ancestry in the controls subset">
##INFO=<ID=non_neuro_AF_amr_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of Latino ancestry in the non_neuro subset">
##INFO=<ID=controls_AF_nfe_nwe,Number=A,Type=Float,Description="Alternate allele frequency in samples of North-Western European ancestry in the controls subset">
##INFO=<ID=controls_AF_nfe_seu,Number=A,Type=Float,Description="Alternate allele frequency in samples of Southern European ancestry in the controls subset">
##INFO=<ID=controls_AF_sas_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of South Asian ancestry in the controls subset">
##INFO=<ID=non_neuro_AF_amr_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of Latino ancestry in the non_neuro subset">
##INFO=<ID=non_cancer_AF_eas_jpn,Number=A,Type=Float,Description="Alternate allele frequency in samples of Japanese ancestry in the non_cancer subset">
##INFO=<ID=non_neuro_AF_nfe_onf,Number=A,Type=Float,Description="Alternate allele frequency in samples of non-Finnish but otherwise indeterminate European ancestry in the non_neuro subset">
##INFO=<ID=non_topmed_AF_eas_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of East Asian ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_afr_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of African-American ancestry in the non_cancer subset">
##INFO=<ID=non_cancer_AF_afr,Number=A,Type=Float,Description="Alternate allele frequency in samples of African-American ancestry in the non_cancer subset">
##INFO=<ID=controls_AF_amr_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of Latino ancestry in the controls subset">
##INFO=<ID=non_neuro_AF_fin_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of Finnish ancestry in the non_neuro subset">
##INFO=<ID=non_neuro_AF_nfe_bgr,Number=A,Type=Float,Description="Alternate allele frequency in samples of Bulgarian ancestry in the non_neuro subset">
##INFO=<ID=non_neuro_AF_oth_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of uncertain ancestry in the non_neuro subset">
##INFO=<ID=non_topmed_AF_nfe_est,Number=A,Type=Float,Description="Alternate allele frequency in samples of Estonian ancestry in the non_topmed subset">
##INFO=<ID=non_topmed_AF_nfe_nwe,Number=A,Type=Float,Description="Alternate allele frequency in samples of North-Western European ancestry in the non_topmed subset">
##INFO=<ID=non_topmed_AF_amr_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of Latino ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_amr,Number=A,Type=Float,Description="Alternate allele frequency in samples of Latino ancestry in the non_cancer subset">
##INFO=<ID=non_topmed_AF_nfe_swe,Number=A,Type=Float,Description="Alternate allele frequency in samples of Swedish ancestry in the non_topmed subset">
##INFO=<ID=non_topmed_AF_nfe_onf,Number=A,Type=Float,Description="Alternate allele frequency in samples of non-Finnish but otherwise indeterminate European ancestry in the non_topmed subset">
##INFO=<ID=controls_AF_eas_kor,Number=A,Type=Float,Description="Alternate allele frequency in samples of Korean ancestry in the controls subset">
##INFO=<ID=non_topmed_AF_eas_oea,Number=A,Type=Float,Description="Alternate allele frequency in samples of non-Korean, non-Japanese East Asian ancestry in the non_topmed subset">
##INFO=<ID=controls_AF_eas_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of East Asian ancestry in the controls subset">
##INFO=<ID=controls_AF_oth_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of uncertain ancestry in the controls subset">
##INFO=<ID=controls_AF_fin,Number=A,Type=Float,Description="Alternate allele frequency in samples of Finnish ancestry in the controls subset">
##INFO=<ID=non_neuro_AF_nfe,Number=A,Type=Float,Description="Alternate allele frequency in samples of non-Finnish European ancestry in the non_neuro subset">
##INFO=<ID=non_neuro_AF_fin_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of Finnish ancestry in the non_neuro subset">
##INFO=<ID=non_cancer_AF_nfe_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of non-Finnish European ancestry in the non_cancer subset">
##INFO=<ID=controls_AF_eas_oea,Number=A,Type=Float,Description="Alternate allele frequency in samples of non-Korean, non-Japanese East Asian ancestry in the controls subset">
##INFO=<ID=non_topmed_AF_nfe_seu,Number=A,Type=Float,Description="Alternate allele frequency in samples of Southern European ancestry in the non_topmed subset">
##INFO=<ID=controls_AF_eas_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of East Asian ancestry in the controls subset">
##INFO=<ID=non_topmed_AF_asj,Number=A,Type=Float,Description="Alternate allele frequency in samples of Ashkenazi Jewish ancestry in the non_topmed subset">
##INFO=<ID=controls_AF_nfe_onf,Number=A,Type=Float,Description="Alternate allele frequency in samples of non-Finnish but otherwise indeterminate European ancestry in the controls subset">
##INFO=<ID=non_topmed_AF_nfe,Number=A,Type=Float,Description="Alternate allele frequency in samples of non-Finnish European ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_oth,Number=A,Type=Float,Description="Alternate allele frequency in samples of uncertain ancestry in the non_cancer subset">
##INFO=<ID=non_topmed_AF_raw,Number=A,Type=Float,Description="Alternate allele frequency in samples in the non_topmed subset, before removing low-confidence genotypes">
##INFO=<ID=non_neuro_AF_nfe_est,Number=A,Type=Float,Description="Alternate allele frequency in samples of Estonian ancestry in the non_neuro subset">
##INFO=<ID=non_topmed_AF_oth_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of uncertain ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_oth_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of uncertain ancestry in the non_cancer subset">
##INFO=<ID=non_cancer_AF_afr_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of African-American ancestry in the non_cancer subset">
##INFO=<ID=non_topmed_AF_afr_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of African-American ancestry in the non_topmed subset">
##INFO=<ID=controls_AF_eas,Number=A,Type=Float,Description="Alternate allele frequency in samples of East Asian ancestry in the controls subset">
##INFO=<ID=non_neuro_AF_eas_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of East Asian ancestry in the non_neuro subset">
##INFO=<ID=non_cancer_AF_nfe_nwe,Number=A,Type=Float,Description="Alternate allele frequency in samples of North-Western European ancestry in the non_cancer subset">
##INFO=<ID=controls_AF_sas,Number=A,Type=Float,Description="Alternate allele frequency in samples of South Asian ancestry in the controls subset">
##INFO=<ID=non_neuro_AF_sas_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of South Asian ancestry in the non_neuro subset">
##INFO=<ID=non_neuro_AF_asj_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of Ashkenazi Jewish ancestry in the non_neuro subset">
##INFO=<ID=non_cancer_AF_nfe_bgr,Number=A,Type=Float,Description="Alternate allele frequency in samples of Bulgarian ancestry in the non_cancer subset">
##INFO=<ID=controls_AF_oth,Number=A,Type=Float,Description="Alternate allele frequency in samples of uncertain ancestry in the controls subset">
##INFO=<ID=non_cancer_AF_eas_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of East Asian ancestry in the non_cancer subset">
##INFO=<ID=non_topmed_AF_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples in the non_topmed subset">
##INFO=<ID=non_neuro_AF_asj,Number=A,Type=Float,Description="Alternate allele frequency in samples of Ashkenazi Jewish ancestry in the non_neuro subset">
##INFO=<ID=non_topmed_AF_eas_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of East Asian ancestry in the non_topmed subset">
##INFO=<ID=non_neuro_AF_raw,Number=A,Type=Float,Description="Alternate allele frequency in samples in the non_neuro subset, before removing low-confidence genotypes">
##INFO=<ID=non_topmed_AF_eas,Number=A,Type=Float,Description="Alternate allele frequency in samples of East Asian ancestry in the non_topmed subset">
##INFO=<ID=non_topmed_AF_fin_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of Finnish ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_asj_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of Ashkenazi Jewish ancestry in the non_cancer subset">
##INFO=<ID=non_topmed_AF_eas_kor,Number=A,Type=Float,Description="Alternate allele frequency in samples of Korean ancestry in the non_topmed subset">
##INFO=<ID=controls_AF_amr_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of Latino ancestry in the controls subset">
##INFO=<ID=non_neuro_AF_eas_oea,Number=A,Type=Float,Description="Alternate allele frequency in samples of non-Korean, non-Japanese East Asian ancestry in the non_neuro subset">
##INFO=<ID=controls_AF_afr_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of African-American ancestry in the controls subset">
##INFO=<ID=controls_AF_amr,Number=A,Type=Float,Description="Alternate allele frequency in samples of Latino ancestry in the controls subset">
##INFO=<ID=non_topmed_AF_eas_jpn,Number=A,Type=Float,Description="Alternate allele frequency in samples of Japanese ancestry in the non_topmed subset">
##INFO=<ID=non_topmed_AF_nfe_bgr,Number=A,Type=Float,Description="Alternate allele frequency in samples of Bulgarian ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_nfe_est,Number=A,Type=Float,Description="Alternate allele frequency in samples of Estonian ancestry in the non_cancer subset">
##INFO=<ID=non_neuro_AF_eas,Number=A,Type=Float,Description="Alternate allele frequency in samples of East Asian ancestry in the non_neuro subset">
##INFO=<ID=non_cancer_AF_nfe,Number=A,Type=Float,Description="Alternate allele frequency in samples of non-Finnish European ancestry in the non_cancer subset">
##INFO=<ID=non_neuro_AF_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples in the non_neuro subset">
##INFO=<ID=non_neuro_AF_sas_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of South Asian ancestry in the non_neuro subset">
##INFO=<ID=controls_AF_nfe_est,Number=A,Type=Float,Description="Alternate allele frequency in samples of Estonian ancestry in the controls subset">
##INFO=<ID=non_topmed_AF_asj_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of Ashkenazi Jewish ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_nfe_swe,Number=A,Type=Float,Description="Alternate allele frequency in samples of Swedish ancestry in the non_cancer subset">
##INFO=<ID=non_topmed_AF_oth,Number=A,Type=Float,Description="Alternate allele frequency in samples of uncertain ancestry in the non_topmed subset">
##INFO=<ID=non_topmed_AF_fin_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of Finnish ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_fin_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of Finnish ancestry in the non_cancer subset">
##INFO=<ID=non_neuro_AF_nfe_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of non-Finnish European ancestry in the non_neuro subset">
##INFO=<ID=controls_AF_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples in the controls subset">
##INFO=<ID=non_cancer_AF_fin,Number=A,Type=Float,Description="Alternate allele frequency in samples of Finnish ancestry in the non_cancer subset">
##INFO=<ID=non_topmed_AF_fin,Number=A,Type=Float,Description="Alternate allele frequency in samples of Finnish ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_eas_oea,Number=A,Type=Float,Description="Alternate allele frequency in samples of non-Korean, non-Japanese East Asian ancestry in the non_cancer subset">
##INFO=<ID=non_topmed_AF_nfe_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of non-Finnish European ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_sas_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of South Asian ancestry in the non_cancer subset">
##INFO=<ID=controls_AF_asj_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of Ashkenazi Jewish ancestry in the controls subset">
##INFO=<ID=non_cancer_AF_raw,Number=A,Type=Float,Description="Alternate allele frequency in samples in the non_cancer subset, before removing low-confidence genotypes">
##INFO=<ID=non_cancer_AF_eas_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of East Asian ancestry in the non_cancer subset">
##INFO=<ID=non_topmed_AF_asj_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of Ashkenazi Jewish ancestry in the non_topmed subset">
##INFO=<ID=non_neuro_AF_oth,Number=A,Type=Float,Description="Alternate allele frequency in samples of uncertain ancestry in the non_neuro subset">
##INFO=<ID=controls_AF_fin_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of Finnish ancestry in the controls subset">
##INFO=<ID=controls_AF_nfe_bgr,Number=A,Type=Float,Description="Alternate allele frequency in samples of Bulgarian ancestry in the controls subset">
##INFO=<ID=controls_AF_asj_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of Ashkenazi Jewish ancestry in the controls subset">
##INFO=<ID=non_topmed_AF_sas_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of South Asian ancestry in the non_topmed subset">
##INFO=<ID=non_cancer_AF_sas,Number=A,Type=Float,Description="Alternate allele frequency in samples of South Asian ancestry in the non_cancer subset">
##INFO=<ID=non_neuro_AF_nfe_seu,Number=A,Type=Float,Description="Alternate allele frequency in samples of Southern European ancestry in the non_neuro subset">
##INFO=<ID=non_cancer_AF_eas_kor,Number=A,Type=Float,Description="Alternate allele frequency in samples of Korean ancestry in the non_cancer subset">
##INFO=<ID=non_topmed_AF_afr_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of African-American ancestry in the non_topmed subset">
##INFO=<ID=controls_AF_sas_male,Number=A,Type=Float,Description="Alternate allele frequency in male samples of South Asian ancestry in the controls subset">
##INFO=<ID=non_topmed_AF_sas_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of South Asian ancestry in the non_topmed subset">
##INFO=<ID=non_topmed_AF_afr,Number=A,Type=Float,Description="Alternate allele frequency in samples of African-American ancestry in the non_topmed subset">
##INFO=<ID=non_neuro_AF_oth_female,Number=A,Type=Float,Description="Alternate allele frequency in female samples of uncertain ancestry in the non_neuro subset">
##INFO=<ID=non_topmed_AF_popmax,Number=A,Type=Float,Description="Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry) in the non_topmed subset">
##INFO=<ID=non_neuro_AF_popmax,Number=A,Type=Float,Description="Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry) in the non_neuro subset">
##INFO=<ID=non_cancer_AF_popmax,Number=A,Type=Float,Description="Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry) in the non_cancer subset">
##INFO=<ID=controls_AF_popmax,Number=A,Type=Float,Description="Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry) in the controls subset">



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
GNOMAD_DIR="/gscmnt/gc2532/dinglab/gnomad_data_2.1/"
#GNOMAD_FILENAME="gnomad.exomes.r2.1.sites.chr"
GNOMAD_FILENAME="gnomad.genomes.r2.1.sites.chr"

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
				#For some mutaations in the VCF there are multiple variants in the variant column. Need to check all fo them against our sites of interest
				variants=data[4].split(',')
				totalvariants=len(variants)
				variantcounter=0
				#initialize lists to store all associated AF numbers
				AF_AFR=[]
				AF_AMR=[]
				AF_ASJ=[]
				AF_EAS=[]
				AF_FIN=[]
				AF_NFE=[]
				AF_OTH=[]
				AF_POPMAX=[]
				for var in variants:
					variantcounter += 1
					vcf_key=data[0]+"."+data[1]+"."+data[3]+"."+var
					if vcf_key in sites_unique:
					 	#If the 2nd variant matches, we need to grab the second AF value in the list, if the 3rd variant matches we need to grab the third time AF matches.
						data2=vline.split(';')
						#Since each VCF line is different we need to go through each element in line and store all the pop allele frequencies based on name.
						for entry in data2:
							if entry.startswith("AF_AFR"):
								AF_AFR=entry.split(',')
							if entry.startswith("AF_AMR"):
								AF_AMR=entry.split(',')
							if entry.startswith("AF_ASJ"):
								AF_ASJ=entry.split(',')
							if entry.startswith("AF_EAS"):
								AF_EAS=entry.split(',')
							if entry.startswith("AF_FIN"):
								AF_FIN=entry.split(',')
							if entry.startswith("AF_NFE"):
								AF_NFE=entry.split(',')
							if entry.startswith("AF_OTH"):
								AF_OTH=entry.split(',')
							if entry.startswith("AF_POPMAX"):
								AF_POPMAX=entry.split(',')
						#	if entry.startswith("AF_FilterStatus"):
						#		AS_FilterStatus=entry.split(',')
							
						actualcounter=variantcounter-1
						print(vcf_key," ",actualcounter,"\n")	
						thing=vcf_key+"\t"+AF_POPMAX[actualcounter]+"\tAF_AFR="+AF_AFR[actualcounter]+"\tAF_AMR="+AF_AMR[actualcounter]+"\tAF_ASJ="+AF_ASJ[actualcounter]+"\tAF_EAS="+AF_EAS[actualcounter]+"\tAF_FIN="+AF_FIN[actualcounter]+"\tAF_NFE="+AF_NFE[actualcounter]+"\tAF_OTH="+AF_OTH[actualcounter]+"\n"
						OUT.write(thing)
					else:
						continue
	fin.close()

OUT.close()

end = time.time()
print(end - start)

####EXAMPLE 2 MUTATIONS IN SAME VCF ENTRY

#zgrep 111957583 /gscmnt/gc2532/dinglab/gnomad_data/gnomad.genomes.r2.0.2.sites.chr1.vcf.bgz
#1	111957583	rs753142160	A	G,ACCCACAGGGGTCAGGGTCTTTTCCCCAGGGCTCACAGACTGATG	5448139.23	PASS	AC=215,1;AF=7.22495e-03,3.36044e-05;AN=29758;BaseQRankSum=-6.45000e-01;ClippingRankSum=-1.16000e-01;DP=673271;FS=1.09400e+00;InbreedingCoeff=7.30000e-03;MQ=6.00000e+01;MQRankSum=-9.20000e-02;QD=1.82800e+01;ReadPosRankSum=-1.81000e-01;SOR=7.81000e-01;VQSLOD=1.55000e+00;VQSR_culprit=SOR;GQ_HIST_ALT=85|1|0|0|0|0|0|0|1|0|1|0|1|0|0|0|1|0|0|213,0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1;DP_HIST_ALT=0|3|13|72|84|65|35|23|4|2|1|0|1|0|0|0|0|0|0|0,0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_HIST_ALT=85|0|1|2|3|2|3|16|28|28|47|34|27|13|8|3|0|0|0|1,0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;GQ_HIST_ALL=408|31|26|30|66|98|95|155|209|295|420|407|965|419|866|480|995|289|837|8403;DP_HIST_ALL=9|35|281|1123|2851|3392|3342|2848|1013|361|131|61|16|18|8|0|0|2|0|0;AB_HIST_ALL=0|9|30|65|202|507|952|1097|1226|741|611|259|121|35|18|5|1|8|39|130;AC_Male=134,1;AC_Female=81,0;AN_Male=16432;AN_Female=13326;AF_Male=8.15482e-03,6.08569e-05;AF_Female=6.07834e-03,0.00000e+00;GC_Male=4457,132,1,1,0,0;GC_Female=3573,79,1,0,0,0;GC_raw=8415,217,2,1,0,0;AC_raw=305,1;AN_raw=30988;GC=8030,211,2,1,0,0;AF_raw=9.84252e-03,3.22706e-05;Hom_AFR=2,0;Hom_AMR=0,0;Hom_ASJ=0,0;Hom_EAS=0,0;Hom_FIN=0,0;Hom_NFE=0,0;Hom_OTH=0,0;Hom=2,0;Hom_raw=2,0;AC_AFR=203,1;AC_AMR=2,0;AC_ASJ=0,0;AC_EAS=2,0;AC_FIN=0,0;AC_NFE=7,0;AC_OTH=1,0;AN_AFR=7932;AN_AMR=824;AN_ASJ=300;AN_EAS=1564;AN_FIN=3460;AN_NFE=14726;AN_OTH=952;AF_AFR=2.55925e-02,1.26072e-04;AF_AMR=2.42718e-03,0.00000e+00;AF_ASJ=0.00000e+00,0.00000e+00;AF_EAS=1.27877e-03,0.00000e+00;AF_FIN=0.00000e+00,0.00000e+00;AF_NFE=4.75350e-04,0.00000e+00;AF_OTH=1.05042e-03,0.00000e+00;STAR_AC=7638;STAR_AC_raw=7880;STAR_Hom=1003;POPMAX=AFR,AFR;AC_POPMAX=203,1;AN_POPMAX=7932,7932;AF_POPMAX=2.55925e-02,1.26072e-04;DP_MEDIAN=23,24;DREF_MEDIAN=3.76188e-41,2.51189e-25;GQ_MEDIAN=99,99;AB_MEDIAN=4.66667e-01,2.08333e-01;AS_RF=8.89651e-01,2.82240e-01;AS_FilterStatus=PASS,RF;CSQ=CCCACAGGGGTCAGGGTCTTTTCCCCAGGGCTCACAGACTGATG|frameshift_variant|HIGH|OVGP1|ENSG00000085465|Transcript|ENST00000369732|protein_coding|11/11||ENST00000369732.3:c.1539_1540insCATCAGTCTGTGAGCCCTGGGGAAAAGACCCTGACCCCTGTGGG|ENSP00000358747.3:p.Tyr514HisfsTer19|1595-1596|1539-1540|513-514|-/HQSVSPGEKTLTPVX|-/CATCAGTCTGTGAGCCCTGGGGAAAAGACCCTGACCCCTGTGGG|rs1126656&COSM4142144&COSM4142145|2||-1||SNV|1|HGNC|8524|YES|||CCDS834.1|ENSP00000358747|Q12889|Q9UJZ3|UPI0000130C53|||||||||||||||AACCCACAGGGGTCAGGGTCTTTTCCCCAGGGCTCACAGACTGATG:0&G:0.001147|AACCCACAGGGGTCAGGGTCTTTTCCCCAGGGCTCACAGACTGATG:1.525e-05&G:6.494e-03|AACCCACAGGGGTCAGGGTCTTTTCCCCAGGGCTCACAGACTGATG:0&G:0.05821|AACCCACAGGGGTCAGGGTCTTTTCCCCAGGGCTCACAGACTGATG:1.553e-05&G:0.006521|AACCCACAGGGGTCAGGGTCTTTTCCCCAGGGCTCACAGACTGATG:0&G:0.0009337|AACCCACAGGGGTCAGGGTCTTTTCCCCAGGGCTCACAGACTGATG:0&G:0|AACCCACAGGGGTCAGGGTCTTTTCCCCAGGGCTCACAGACTGATG:0&G:0.0006633|AACCCACAGGGGTCAGGGTCTTTTCCCCAGGGCTCACAGACTGATG:0&G:0.008333||0&1&1|0&1&1||||||HC|||POSITION:0.755522827687776&PHYLOCSF_TOO_SHORT
