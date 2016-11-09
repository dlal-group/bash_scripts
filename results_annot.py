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

#include a usage message
parser=argparse.ArgumentParser()
parser.add_argument('<result file>')
parser.add_argument('<annotation_vcf>')
parser.add_argument('<GEMMA/genABEL>')
if len(sys.argv)==1:
    parser.print_usage()
    sys.exit(1)
args=parser.parse_args()

res_file=sys.argv[1]
annot_file=sys.argv[2]
mode=sys.argv[3]


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
sys.stderr.write('Annotation read in '+ elapsed_time+'...\n')

sys.stderr.write('Starting annotation...\n')
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
	print 'SNP,Chromosome,Position,A0,A1,NoMeasured,CallRate,Pexact,MarkerType,Rsq,p,beta,sebeta,effallelefreq,MAF,strand,rsID'
	with open('%s' %(res_file) ,'r') as current_file:
		next(current_file)
		for line in current_file:
			site=line.rstrip().split(",")
			site_key="_".join([site[1],site[2],site[3],site[4]])
			try:
				all_annots[site_key]
				# all_res[site_key]=[site,all_annots[site_key]]
				print '%s,%s' %(",".join(site), all_annots[site_key])
			except KeyError, e:
				# all_res[site_key]=[site,":".join([site[1],site[2]])]
				print '%s,%s' %(",".join(site),":".join([site[1],site[2]]))
