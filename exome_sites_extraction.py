#!/usr/bin/env python2.7
#Script used to extract exome chip sites from:
# 1) exome chip
# 2) 1000G data
# 3) imputed data

#  and calculate frequencies with allele consistency match/check

import gzip 
import re 
import sys 
import subprocess as sub
import time
import collections
import itertools
from itertools import chain
import os

#args needed:
# - single chrom exome chip bim file
# - single chrom vcf file from 1000G
# - single chrom info file from imputation

exome_chip=sys.argv[1]
tgp_data=sys.argv[2]
impute_data=sys.argv[3]
chromosome=sys.argv[4]
pop=sys.argv[5]
outpath=sys.argv[6]

# exome_chip="/netapp/dati/WGS_REF_PANEL/genotypes/CARL/exome/2/chr2_flipped.bim"
# tgp_data="/home/cocca/imputation/31012017_MERGED_TEST/CARL/2_exome_TGP_sites.txt"
# impute_data="/home/cocca/imputation/MERGED_INGI_TGP3_23012017/CARL/MERGED/ALL/chr2.gen_info"
# chromosome=2
# pop="CARL"
# outpath="/home/cocca/imputation/31012017_MERGED_TEST/CARL"
#set constant values
carl_samples=630
fvg_samples=1590

if pop == "CARL":
	current_samples=carl_samples
elif pop == "FVG":
	current_samples=fvg_samples

#first we need to read the exome chip bim file
exome_data={}
with open('%s' %(exome_chip) ,'r') as current_exome:
	for line in current_exome:
		site=line.rstrip().split("\t")
		chrom=site[0]
		pos=site[3]
		rsID=site[1]
		a0=site[4]
		a1=site[5]
		site_k=pos
		exome_data[site_k]= {"chrom":chrom,"pos":pos,"rsID":rsID,"a0":a0,"a1":a1}
		#generate a list to be used with bcftools to extract data by position, 
		# so we don't have to worry about the vcf format file


tgp_values={}
#read the data from TGP: we'll get the data from the table already created
with open('%s' %(tgp_data) ,'r') as current_tgp:
	for line in current_tgp:
		site=line.rstrip().split("\t")
		chrom=site[0]
		pos=site[1]
		rsID=site[4]
		a0=site[2]
		a1=site[3]
		a1_af=site[7]
		an=int(site[5])
		ac1=int(site[6])
		ac0=an - ac1
		a0_af=1-float(a1_af)
		site_k=chrom+"_"+pos+"_"+a0+"_"+a1
		tgp_values[site_k]= {"chrom":chrom,"pos":pos,"rsID":rsID,"a0":a0,"a1":a1,"a1_af":a1_af,"a0_af":a0_af,"an":an,"ac1":ac1,"ac0":ac0}


#read the info file
info_values={}
with open('%s' %(impute_data) ,'r') as current_info:
	next(current_info)
	for line in current_info:
		site=line.rstrip().split(" ")
		chrom=chromosome
		pos=site[2]
		rsID=site[1]
		a0=site[3]
		a1=site[4]
		a1_af=float(site[5])
		a0_af=1-a1_af
		an=current_samples*2
		ac1=a1_af*an
		ac0=a0_af*an
		site_k=str(chrom)+"_"+pos+"_"+a0+"_"+a1
		info_values[site_k]= {"chrom":chrom,"pos":pos,"rsID":rsID,"a0":a0,"a1":a1,"a1_af":a1_af,"a0_af":a0_af,"an":an,"ac1":ac1,"ac0":ac0}


#now try to get from the exome chip data, the imputed one
#chek for allele flipping checking for the key with a0_a1 and if not with a1_a0
all_imputed_stuff={}
flipped_imputed_stuff={}
missing_imputed_stuff={}
for exome_pos in exome_data:
	exome_site=exome_data[exome_pos]
	check_key1="_".join([exome_site["chrom"],exome_site["pos"],exome_site["a0"],exome_site["a1"]]) 
	check_key2="_".join([exome_site["chrom"],exome_site["pos"],exome_site["a1"],exome_site["a0"]]) 
	try:
		info_values[check_key1]
		# all_res[site_key]=[site,all_annots[site_key]]
	except KeyError, e:
		# all_res[site_key]=[site,":".join([site[1],site[2]])]
		try:
			info_values[check_key2]
		except KeyError, e:
			missing_imputed_stuff[check_key1]= exome_data[exome_pos]
		else:
			flipped_imputed_stuff[check_key2] = info_values[check_key2]
			all_imputed_stuff[check_key2] = info_values[check_key2]
	else:
		all_imputed_stuff[check_key1] = info_values[check_key1]

all_data_to_check={}
missing_to_check={}
flipped_to_check={}
final_compare={}
#now get the data from TGP
for check_site in all_imputed_stuff:
	t_check_key1=check_site
	t_check_key2="_".join([check_site.split("_")[0],check_site.split("_")[1],check_site.split("_")[3],check_site.split("_")[2]])
	try:
		tgp_values[t_check_key1]
		# all_res[site_key]=[site,all_annots[site_key]]
	except KeyError, e:
		# all_res[site_key]=[site,":".join([site[1],site[2]])]
		try:
			tgp_values[t_check_key2]
		except KeyError, e:
			missing_to_check[t_check_key1]= all_imputed_stuff[check_site]
		else:
			flipped_to_check[t_check_key2] = tgp_values[t_check_key2]
			all_data_to_check[t_check_key2] = tgp_values[t_check_key2]
	else:
		all_data_to_check[t_check_key1] = tgp_values[t_check_key1]
		final_compare[t_check_key1] = {"chrom":tgp_values[t_check_key1]["chrom"],"pos":tgp_values[t_check_key1]["pos"],"rsID":tgp_values[t_check_key1]["rsID"],"a0":all_imputed_stuff[t_check_key1]["a0"],"a1":all_imputed_stuff[t_check_key1]["a1"],"tgp_a1_f":tgp_values[t_check_key1]["chrom"],"tgp_a1_f":tgp_values[t_check_key1]["a1_af"],"tgp_an":tgp_values[t_check_key1]["an"],"tgp_ac1":tgp_values[t_check_key1]["ac1"],"imp_a1_f":all_imputed_stuff[t_check_key1]["a1_af"],"imp_an":all_imputed_stuff[t_check_key1]["an"],"imp_ac1":all_imputed_stuff[t_check_key1]["ac1"]}

#write a file for testing allele frequencies
freq_test=open('%s/chr%s_to_test.txt' %(outpath, chrom),'w')
print >> freq_test,"chrom\tpos\trsID\ta0\ta1\ttgp_an\ttgp_ac1\ttgp_a1_f\timp_an\timp_ac1\timp_a1_f"

for to_check_site in final_compare:
	site_to_check=final_compare[to_check_site]
	print >> freq_test,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(site_to_check['chrom'], site_to_check['pos'],site_to_check['rsID'], site_to_check['a0'], site_to_check['a1'], site_to_check['tgp_an'], site_to_check['tgp_ac1'], site_to_check['tgp_a1_f'], site_to_check['imp_an'], site_to_check['imp_ac1'], site_to_check['imp_a1_f'])