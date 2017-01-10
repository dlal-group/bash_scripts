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
import select
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool 


# if select.select([sys.stdin,],[],[],0.0)[0]:
# 	res_file=sys.stdin # res_file="/home/cocca/analyses/1000G_test/CARL/MCH_out/CARL_MCH_10_Oct_08_2016_cocca_results.csv"
# 	res_file="/home/cocca/imputation/MERGED_DATA/CARL/MERGED/ALL/chr20.gen_info"
# 	annot_file=sys.argv[1] # annot_file="/netapp/nfs/resources/1000GP_phase3/vcf/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
# 	mode=sys.argv[2] # mode="genABEL"
# 	#this argument is optional and used only in genABEL mode
# 	i_conv=sys.argv[3] # i_conv="/netapp02/data/imputation/INGI_TGP3/CARL/carl/MERGED/CLEANED/chr10.gen_info"
# 	sys.stderr.write('annot file: '+annot_file+'\n')
# else:
#include a usage message
parser=argparse.ArgumentParser()
parser.add_argument('<result file>')
parser.add_argument('<annotation_vcf>')
parser.add_argument('<GEMMA/genABEL/datABEL/info>')
parser.add_argument('[indel recode file]', nargs='?', default=None)
# parser.add_argument('[threads]')
if len(sys.argv)==1:
    parser.print_usage()
    sys.exit(1)
args=parser.parse_args()
res_file=sys.argv[1] # res_file="/home/cocca/analyses/1000G_test/CARL/MCH_out/CARL_MCH_10_Oct_08_2016_cocca_results.csv"
# res_file="/netapp02/data/imputation/INGI_TGP3/VBI/dose/VBI_INGI_TGP3_chr10.map"
annot_file=sys.argv[2] # annot_file="/netapp/nfs/resources/1000GP_phase3/vcf/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.180000.vcf.gz"
mode=sys.argv[3] # mode="datABEL"
#this argument is optional and used only in genABEL mode
# i_conv=sys.argv[4] # i_conv="/netapp02/data/imputation/INGI_TGP3/CARL/carl/MERGED/CLEANED/chr10.gen_info"
# threads=sys.argv[5]#threads=8

def ann_reader(d,ann_line):
	x=ann_line.split("\t")
	if not re.match('#',line):
		ann_key="_".join([x[0],x[1],x[3],x[4]])
		d[ann_key] = x[2]

	# return d

def recode_reader(info_file):
	recoded={}
	with open('%s' %(info_file) ,'r') as current_recode:
		next(current_recode)
		# if not re.match('SNP',line):
		for line in current_recode:
			if line.rstrip().split(" ")[0] == "---":
				rs_id=line.rstrip().split(" ")[1]
				rec_key="_".join((rs_id.split("_")[0].split(":")))
				recoded[rec_key] = "_".join(rs_id.split(":"))
			else:
				rec_key="_".join((line.rstrip().split(" ")[0],line.rstrip().split(" ")[2]))
				recoded[rec_key] = "_".join((rec_key,line.rstrip().split(" ")[3],line.rstrip().split(" ")[4]))
	return recoded


# read annotation file (it will be in vcf.gz format)
sys.stderr.write('reading annotation file...\n')
start_time = time.time()

# pool=ThreadPool(int(threads))
current_annot=gzip.open('%s' %(annot_file), 'r')
all_annots={}

for line in current_annot:
	x=line.split()
	if not re.match('#',line):
		ann_key="_".join([x[0],x[1],x[3],x[4]])
		all_annots[ann_key] = x[2]

# for line in current_annot:
# 	pool.map(ann_reader(all_annots,line),line)
elapsed_time = time.time() - start_time
sys.stderr.write('Annotation read in '+ str(timedelta(seconds=elapsed_time)) +'...\n')

sys.stderr.write('Starting annotation...\n')
start_time = time.time()
#now we need to get from res files all keys and add the rsID
all_res = {}
#read res_file
if mode == 'GEMMA':
	print 'chr\trs\tps\tn_mis\tn_obs\tallele1\tallele0\taf\tbeta\tse\tp_wald\tp_lrt\tp_score\trsID'
	if select.select([sys.stdin,],[],[],0.0)[0]:
		current_file=sys.stdin
	else:
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
	all_recoded=recode_reader(i_conv)
	
	# all_recoded={}
	# with open('%s' %(i_conv) ,'r') as current_recode:
	# 	next(current_recode)
	# 	# if not re.match('SNP',line):
	# 	for line in current_recode:
	# 		if line.rstrip().split(" ")[0] == "---":
	# 			rs_id=line.rstrip().split(" ")[1]
	# 			rec_key="_".join((rs_id.split("_")[0].split(":")))
	# 			all_recoded[rec_key] = "_".join(rs_id.split(":"))
	# 		else:
	# 			rec_key="_".join((line.rstrip().split(" ")[0],line.rstrip().split(" ")[2]))
	# 			all_recoded[rec_key] = "_".join((rec_key,line.rstrip().split(" ")[3],line.rstrip().split(" ")[4]))

	elapsed_time_conv = time.time() - start_time_conv
	sys.stderr.write('Conversion file read in '+ str(timedelta(seconds=elapsed_time_conv)) +'...\n')

	print 'SNP,Chromosome,Position,A0,A1,NoMeasured,CallRate,Pexact,MarkerType,Rsq,p,beta,sebeta,effallelefreq,MAF,strand,rsID'
	with open('%s' %(res_file) ,'r') as current_file:
		for line in current_file:
			if not re.match('SNP',line):
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

elif mode == 'datABEL':
	#we need to annotate imputed files in map/info format
	sys.stderr.write('Reading conversion file for datABEL format...\n')
	start_time_conv = time.time()
	# we need to retrieve the indel recoded conversion
	all_recoded=recode_reader(i_conv)
	
	elapsed_time_conv = time.time() - start_time_conv
	sys.stderr.write('Conversion file read in '+ str(timedelta(seconds=elapsed_time_conv)) +'...\n')

	print 'SNP Position A0 A1 Rsq'
	with open('%s' %(res_file) ,'r') as current_file:
		for line in current_file:
			if not re.match('SNP',line):
				site=line.rstrip().split(" ")
				site_key="_".join([site[0].split(":")[0],site[1]])
				# site_key="_".join([site[1],site[2],site[3],site[4]])
				try:
					all_annots[all_recoded[site_key]]
					# all_res[site_key]=[site,all_annots[site_key]]
					print '%s %s %s %s %s' %(all_annots[all_recoded[site_key]],site[1],site[2],site[3],site[4])
				except KeyError, e:
					# all_res[site_key]=[site,":".join([site[1],site[2]])]
					print '%s' %(" ".join(site))

elif mode == 'info':
	print 'snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0'
	with open('%s' %(res_file) ,'r') as current_file:
			for line in current_file:
				if not re.match('snp_id',line):
					site=line.rstrip().split(" ")
					if re.match('---',site[0]):
						site_key="_".join([site[1].split(":")[0],site[1].split(":")[1]])
					else:
						site_key="_".join([site[0],site[2],site[3],site[4]])
						
					try:
						all_annots[site_key]
						# all_res[site_key]=[site,all_annots[site_key]]
						print '%s %s %s %s %s %s %s %s %s %s %s %s' %(site[0],all_annots[site_key],site[2],site[3],site[4],site[5],site[6],site[7],site[8],site[9],site[10],site[11])
					except KeyError, e:
						# all_res[site_key]=[site,":".join([site[1],site[2]])]
						print '%s' %(" ".join(site))


elapsed_time = time.time() - start_time

sys.stderr.write('Annotation done in '+ str(timedelta(seconds=elapsed_time)) +'...\n')
