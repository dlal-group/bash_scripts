#!/usr/bin/env python2.7
#annotate variants after genabel analyses using 1000GPh3 data
import argparse
import gzip 
import re 
import sys 
import subprocess as sub
import time
import collections
import itertools
from itertools import chain
from datetime import timedelta

#include a usage message
parser=argparse.ArgumentParser()
parser.add_argument('<result file>')
parser.add_argument('<annotation_vcf>')
parser.add_argument('<GEMMA/genABEL>')
parser.add_argument('[indel recode file]')
if len(sys.argv)==1:
    parser.print_usage()
    sys.exit(1)
args=parser.parse_args()

res_file=sys.argv[1] # res_file="/home/cocca/analyses/1000G_test/CARL/MCH_out/CARL_MCH_10_Oct_08_2016_cocca_results.csv"
annot_file=sys.argv[2] # annot_file="/netapp/nfs/resources/1000GP_phase3/vcf/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.180000.vcf.gz"
mode=sys.argv[3] # mode="genABEL"
#this argument is optional and used only in genABEL mode
i_conv=sys.argv[4] # i_conv="/netapp02/data/imputation/INGI_TGP3/CARL/carl/MERGED/CLEANED/chr10.gen_info"


# read annotation file (it will be in vcf.gz format)
current_annot=gzip.open('%s' %(annot_file), 'r')
sys.stderr.write('reading annotation file...\n')
start_time = time.time()

all_annots={}
for line in current_annot:
	x=line.split()
	if not re.match('#',line):
		ann_key="_".join([x[0],x[1],x[3],x[4]])
		all_annots[ann_key] = x[2]

elapsed_time = time.time() - start_time

sys.stderr.write('Annotation read in '+ str(timedelta(seconds=elapsed_time)) +'...\n')

sys.stderr.write('Starting annotation...\n')
start_time = time.time()
#now we need to get from res files all keys and add the rsID
all_res = {}
#read res_file
if mode == 'GEMMA':
	print 'chr\trs\tps\tn_mis\tn_obs\tallele1\tallele0\taf\tbeta\tse\tp_wald\tp_lrt\tp_score\trsID'
	current_file=gzip.open('%s' %(res_file), 'r')
	for line in current_file:
		if not re.match('chr',line):		
			site=line.rstrip().split("\t")
			site_key="_".join([site[0],site[2],site[5],site[6]])
			try:
				all_annots[site_key]
				# all_res[site_key]=[site,all_annots[site_key]]
				print '%s\t%s' %("\t".join(site), all_annots[site_key])
			except KeyError, e:
				# all_res[site_key]=[site,":".join([site[0],site[2]])]
				print '%s\t%s' %("\t".join(site),":".join([site[0],site[2]]))

elif mode == 'genABEL':
	#TODO: add ability to read result csv from GWA output from tar file with
	# tar -axf VB_TG_22_Oct_08_2016_cocca.tar.gz VB_TG_22_Oct_08_2016_cocca.log --to-stdout
	sys.stderr.write('Reading conversion file for genABEL format...\n')
	start_time_conv = time.time()
	# we need to retrieve the indel recoded conversion
	all_recoded={}
	with open('%s' %(i_conv) ,'r') as current_recode:
		next(current_recode)
		for line in current_recode:
			if line.rstrip().split(" ")[0] == "---":
				rs_id=line.rstrip().split(" ")[1]
				rec_key="_".join((rs_id.split("_")[0].split(":")))
				all_recoded[rec_key] = "_".join(rs_id.split(":"))
			else:
				rec_key="_".join((line.rstrip().split(" ")[0],line.rstrip().split(" ")[2]))
				all_recoded[rec_key] = "_".join((rec_key,line.rstrip().split(" ")[3],line.rstrip().split(" ")[4]))

	elapsed_time_conv = time.time() - start_time_conv
	sys.stderr.write('Conversion file read in '+ str(timedelta(seconds=elapsed_time_conv)) +'...\n')

	print 'SNP,Chromosome,Position,A0,A1,NoMeasured,CallRate,Pexact,MarkerType,Rsq,p,beta,sebeta,effallelefreq,MAF,strand,rsID'
	with open('%s' %(res_file) ,'r') as current_file:
		next(current_file)
		for line in current_file:
			site=line.rstrip().split(",")
			site_key="_".join([site[1],site[2]])
			# site_key="_".join([site[1],site[2],site[3],site[4]])
			try:
				all_annots[all_recoded[site_key]]
				# all_res[site_key]=[site,all_annots[site_key]]
				print '%s,%s' %(",".join(site), all_annots[all_recoded[site_key]])
			except KeyError, e:
				# all_res[site_key]=[site,":".join([site[1],site[2]])]
				print '%s,%s' %(",".join(site),":".join([site[1],site[2]]))

elapsed_time = time.time() - start_time

sys.stderr.write('Annotation done in '+ str(timedelta(seconds=elapsed_time)) +'...\n')
