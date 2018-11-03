#Reyka Jayasinghe (reyka@wustl.edu)
#Edited: October 16th, 2018
#Usage: python gnomad_annotation.py mutationfile
#Note: Only works for snps, if indels are in mutation file it will throw an error

#To download gnomad data
#http://gnomad.broadinstitute.org/downloads
#mkdir gnomad_data

#Download relevant data: gsutil cp gs://gnomad-public/release/2.1/vcf/exomes/* .

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
#Site	controls_AF_afr_male	non_neuro_AF_eas_kor	non_topmed_AF_amr	non_cancer_AF_asj_female	non_cancer_AF_oth_female	non_neuro_AF_asj_female	non_neuro_AF_afr_male	controls_AF_nfe_swe	non_neuro_AF_afr_femalenon_topmed_AF_amr_female	non_cancer_AF_female	non_cancer_AF_nfe_onf	non_cancer_AF_male	non_topmed_AF_oth_female	non_cancer_AF_sas_female	non_neuro_AF_female	controls_AF_afr	non_neuro_AF_eas_jpn	non_cancer_AF_amr_male	controls_AF_fin_male	non_neuro_AF_nfe_nwe	non_topmed_AF_nfe_male	non_neuro_AF_sas	non_cancer_AF_fin_male	non_cancer_AF_nfe_seu	non_neuro_AF_nfe_female	non_neuro_AF_afr	controls_AF_raw	non_cancer_AF_eas	non_cancer_AF_amr_female	non_neuro_AF_nfe_swe	controls_AF_male	non_topmed_AF_male	controls_AF_eas_jpn	controls_AF_nfe_female	non_neuro_AF_amr	non_neuro_AF_eas_female	controls_AF_nfe_male	non_neuro_AF_fin	non_topmed_AF_sas	non_cancer_AF_nfe_female	non_cancer_AF_asj	controls_AF_nfe	controls_AF_oth_female	controls_AF_asj	non_neuro_AF_amr_male	controls_AF_nfe_nwe	controls_AF_nfe_seu	controls_AF_sas_female	non_neuro_AF_amr_female	non_cancer_AF_eas_jpn	non_neuro_AF_nfe_onf	non_topmed_AF_eas_male	non_cancer_AF_afr_male	non_cancer_AF_afr	controls_AF_amr_female	non_neuro_AF_fin_male	non_neuro_AF_nfe_bgr	non_neuro_AF_oth_male	non_topmed_AF_nfe_est	non_topmed_AF_nfe_nwe	non_topmed_AF_amr_male	non_cancer_AF_amr	non_topmed_AF_nfe_swe	non_topmed_AF_nfe_onf	controls_AF_eas_kor	non_topmed_AF_eas_oea	controls_AF_eas_male	controls_AF_oth_male	controls_AF_fin	non_neuro_AF_nfe	non_neuro_AF_fin_female	non_cancer_AF_nfe_male	controls_AF_eas_oea	non_topmed_AF_nfe_seu	controls_AF_eas_female	non_topmed_AF_asj	controls_AF_nfe_onf	non_topmed_AF_nfe	non_cancer_AF_oth	non_topmed_AF_raw	non_neuro_AF_nfe_est	non_topmed_AF_oth_male	non_cancer_AF_oth_male	non_cancer_AF_afr_female	non_topmed_AF_afr_male	controls_AF_eas	non_neuro_AF_eas_male	non_cancer_AF_nfe_nwe	controls_AF_sas	non_neuro_AF_sas_male	non_neuro_AF_asj_male	non_cancer_AF_nfe_bgr	controls_AF_oth	non_cancer_AF_eas_female	non_topmed_AF_female	non_neuro_AF_asj	non_topmed_AF_eas_female	non_neuro_AF_raw	non_topmed_AF_eas	non_topmed_AF_fin_male	non_cancer_AF_asj_male	non_topmed_AF_eas_kor	controls_AF_amr_male	non_neuro_AF_eas_oea	controls_AF_afr_female	controls_AF_amr	non_topmed_AF_eas_jpn	non_topmed_AF_nfe_bgr	non_cancer_AF_nfe_est	non_neuro_AF_eas	non_cancer_AF_nfe	non_neuro_AF_male	non_neuro_AF_sas_female	controls_AF_nfe_est	non_topmed_AF_asj_female	non_cancer_AF_nfe_swe	non_topmed_AF_oth	non_topmed_AF_fin_female	non_cancer_AF_fin_female	non_neuro_AF_nfe_male	controls_AF_female	non_cancer_AF_fin	non_topmed_AF_fin	non_cancer_AF_eas_oea	non_topmed_AF_nfe_female	non_cancer_AF_sas_male	controls_AF_asj_male	non_cancer_AF_raw	non_cancer_AF_eas_male	non_topmed_AF_asj_male	non_neuro_AF_oth	controls_AF_fin_female	controls_AF_nfe_bgr	controls_AF_asj_female	non_topmed_AF_sas_male	non_cancer_AF_sas	non_neuro_AF_nfe_seu	non_cancer_AF_eas_kor	non_topmed_AF_afr_female	controls_AF_sas_male	non_topmed_AF_sas_female	non_topmed_AF_afr	non_neuro_AF_oth_female	non_topmed_AF_popmax	non_neuro_AF_popmax	non_cancer_AF_popmax	controls_AF_popmax
#1.120298095.A.C	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	1.59084e-05	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	7.48817e-06	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	3.59583e-05	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	3.31170e-05	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	1.11662e-05	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	8.95480e-06	0.00000e+00	4.08367e-06	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	4.80455e-06	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	8.89142e-06	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	2.01110e-05	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	8.95480e-06	1.11662e-05	NA	NA


#Import necessary packages
from collections import defaultdict
import sys
import gzip
import time

start = time.time()

#Directory where Gnomad Data is located
#File name scheme: gnomad.genomes.r2.0.2.sites.chr1.vcf.bgz
GNOMAD_DIR="/gscmnt/gc2532/dinglab/gnomad_data_2.1/"
GNOMAD_FILENAME="gnomad.exomes.r2.1.sites.chr"

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

#Dictionary to store values
complicatedAF=defaultdict(dict)
AF_types=[]
all_mutations=[]

for line in MUT:
	(chromosome,position,refallele,altallele) = line.strip().split('\t')
#Skip all indels
	if len(refallele) > 1:
		continue
	if len(altallele) > 1:
		continue
	key = chromosome+"."+position+"."+refallele+"."+altallele
	if key not in all_mutations:
		all_mutations.append(key)
	mutations[chromosome].append(key)
		
MUT.close()

for chrom,sites in mutations.items():
	processed=[]
	#Get unique set of sites
	sites_unique=set(sites)
	#open associated gnomad VCF file to find allele frequency information
	VCF=GNOMAD_DIR+GNOMAD_FILENAME+chrom+".vcf.bgz"
	print(chrom,":",sites,",",VCF)
	with gzip.open(VCF,'r') as fin:
		for vline in fin:
			if not vline.startswith("#"):
				data=vline.strip('\'').split('\t')
				#For some mutaations in the VCF there are multiple variants in the variant column. Need to check all of them against our sites of interest
				variants=data[4].split(',')
				totalvariants=len(variants)
				variantcounter=0
				for var in variants:
					#all variants for this cancer type have been processed
					if set(processed) == set(sites_unique):
						break
					else:
						pass
					variantcounter += 1
					vcf_key=data[0]+"."+data[1]+"."+data[3]+"."+var
					if vcf_key in sites_unique:
					 	#If the 2nd variant matches, we need to grab the second AF value in the list, if the 3rd variant matches we need to grab the third time AF matches.
						data2=vline.split(';')
						#Since each VCF line is different we need to go through each element in line and store all the pop allele frequencies based on name.
						for entry in data2:
							if "_AF_" in entry:
								actualcounter=variantcounter-1
								(type_AF,AFvalues)=entry.split('=')
								values=AFvalues.split(',')							
								complicatedAF[vcf_key][type_AF]=values[actualcounter]
								if type_AF not in AF_types:
									AF_types.append(type_AF)
						#check if all variants have been processed for this cancer type
						processed.append(vcf_key)
					else:
						continue
			if set(processed) == set(sites_unique):
				break	
			else:
				pass
	fin.close()

#Print header
header="Site"+"\t"+"\t".join(AF_types)+"\n"
OUT.write(header)


#########Print all values
for mut in all_mutations:
	print("Processing:",mut)
	mutationdata=mut
	OUT.write(mutationdata)
	for af_val in AF_types:
		try:
			freq=complicatedAF[mut][af_val]
		except KeyError:
			freq="NA"
		AFdata="\t"+freq
		OUT.write(AFdata)
	OUT.write("\n")

OUT.close()

end = time.time()
print(end - start)
