#!/usr/bin/env python2.7
#script used to re filter vcf files for reference generation
import gzip 
import re 
import sys 
import subprocess as sub
import time
import collections
import itertools


#prepare data in tab format to be used to add annotations in VCF files
#right now we have six colunms TODO: extend to different number of cols

annotdata = sys.stdin.readlines()
# annotdata ="test_annots_dupe.tsv"

# Easy peasy:create a dictionary with position and ref alelle as key
# Values will be annotations

all_sites=collections.defaultdict(lambda: collections.defaultdict(list))
with open('%s' %(annotdata) ,'r') as anno_file:
	next(anno_file)
	for c_row in anno_file:
		site=c_row.rstrip().split(" ")
		if len(site[2]) != len(site[3]):
			all_sites[site[1]]["INDEL"].append([site[0],site[1],site[2],site[3],site[4],site[5]])
		else:
			all_sites[site[1]]["SNP"].append([site[0],site[1],site[2],site[3],site[4],site[5]])

#We need to collapse overlapping annotations

for key in all_sites:
# key='119902046'
	for v_type in all_sites[key]:
		#we need to behave in a different way if it's multiallelic
		if len(all_sites[key][v_type]) == 1:
			variant = map(str,all_sites[key][v_type][0])
			print "\t".join(variant)
			# print variant[0],' ',variant[1],' ',variant[2],' ',variant[3],' ',v_type
		else:
			#we need to collapse stuff and remove duplicates and print
			collapsed = zip(*all_sites[key][v_type])
			collapsed_all=[]
			for a in collapsed:
				collapsed_all.append(list(set(a)))
			print "\t".join([",".join(sublist) for sublist in collapsed_all])