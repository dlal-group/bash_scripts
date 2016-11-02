#!/usr/bin/env python2.7
#Script used to merge data in a isec style
import gzip 
import re 
import sys 
import subprocess as sub
import time
import collections
import itertools
from itertools import chain

all_files=sys.argv[1:3]
# file2=sys.argv[2]
# file3=sys.argv[3]
pops=sys.argv[4:]

# pops=["CARL","FVG","VBI"]
# file1="/netapp/dati/INGI_WGS/18112015/"+pops[0]+"/12112015_FILTERED_REL/30092016_UNRELATED/ALL_"+pops[0]+"_02102016.vcf.gz.freq.tab.10000"
# file2="/netapp/dati/INGI_WGS/18112015/"+pops[1]+"/12112015_FILTERED_REL/30092016_UNRELATED/ALL_"+pops[1]+"_02102016.vcf.gz.freq.tab.10000"
# file3="/netapp/dati/INGI_WGS/18112015/"+pops[2]+"/12112015_FILTERED_REL/30092016_UNRELATED/ALL_"+pops[2]+"_02102016.vcf.gz.freq.tab.10000"

# all_files=[file1,file2,file3]
# all_sites_1=collections.defaultdict(lambda: collections.defaultdict(list))
# all_sites_2=collections.defaultdict(lambda: collections.defaultdict(list))
# all_sites_3=collections.defaultdict(lambda: collections.defaultdict(list))
pop_files={}
for file in all_files:
	pop = all_files.index(file)
	pop_files[pops[pop]] = file

all_sites=collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))
for k in pop_files:
	#read each file and create a dictionary with chr_pos_ref_alt as key and all other fields as values
	with open('%s' %(pop_files[k]) ,'r') as current_file:
		# for line in gzip.open(file1, 'r'):
		next(current_file)
		for line in current_file:
			site=line.rstrip().split("\t")
			site_k="_".join([site[0],site[1]])
			if len(site[3]) != len(site[4]):
				all_sites[site_k]["INDEL"][k].append([site[2:]])
			else:
				all_sites[site_k]["SNP"][k].append([site[2:]])

#now we need to check the overlap between each population
# 1_1000894
print 'CHROM\tPOS\tID\tREF\tALT',
for pop in pops:
	print '\t%s\t%s\t%s\t%s' %(pop+"_AC",pop+"_AN",pop+"_AAF",pop+"_MAF"),
print '\tPOPS'

for var in all_sites:
	for v_type in all_sites[var]:
		f_k=[all_sites[var][v_type].keys()][0][0]
		var_split=var.split("_")
		var_info=list(chain.from_iterable(all_sites[var][v_type][f_k][0]))
		print '%s\t%s\t%s\t%s\t%s' %(var_split[0],var_split[1],var_info[0],var_info[1],var_info[2]),
		for pop in pops:
			if pop in all_sites[var][v_type].keys():
				info=list(chain.from_iterable(all_sites[var][v_type][pop][0]))
				print '\t%s\t%s\t%s\t%s' %(info[3],info[4],info[5],info[6]),
			else:
				print '\tNA\tNA\tNA\tNA',
		print '\t%s' %("_".join(all_sites[var][v_type].keys()))
		# print '\r'