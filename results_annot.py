#!/usr/bin/env python2.7
#annotate variants after genabel analyses using 1000GPh3 data
import gzip 
import re 
import sys 
import subprocess as sub
import time
import collections
import itertools
from itertools import chain



res_file=sys.argv[1]
annot_file=sys.argv[2]


# read annotation file (it will be in vcf.gz format)
current_annot=gzip.open('%s' %(annot_file), 'r')

all_annots={}
for line in current_annot:
	x=line.split()
	if not re.match('#',line):
		ann_key="_".join([x[0],x[1],x[3],x[4]])
		all_annots[ann_key] = x[2]

#now we need to get from res files all keys and add the rsID
all_res = {}
#read res_file
with open('%s' %(res_file) ,'r') as current_file:
	next(current_file)
	for line in current_file:
		site=line.rstrip().split(",")
		site_key="_".join([site[1],site[2],site[3],site[4]])
		if site_key in all_annots.keys():
			all_res[site_key]=[site,all_annots[site_key]]
			print '%s,%s' %(",".join(site), all_annots[site_key])
		else:
			all_res[site_key]=[site,":".join([site[1],site[2]])]
			print '%s,%s' %(",".join(site),":".join([site[1],site[2]]))



