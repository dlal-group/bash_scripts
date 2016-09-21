#!/usr/bin/env python2.7
#script used to re filter vcf files for reference generation
import re 
import sys 
import subprocess as sub
import time
import collections
import random

#################################################
# Script to randomly select a certain percentage of sites to remove from a list based on MAF
#INPUT -> frequency tab calculated by plink and reformatted by awk
#OUTPUT -> list of sites to remove

start_time = time.time()
start_time1 = time.ctime(int(start_time))
print start_time1

#freq table
freq_tab=sys.argv[1]
# freq_tab="chr2_freq.tab"

# out_list=sys.argv[2]
# out_list="chr2_freq.tab.exclude"

#define frequecy bins with upper and lower limits
maf_bins={1:[0,0.1], 2:[0.1,0.2],3:[0.2,0.3],4:[0.3,0.4],5:[0.4,0.5],6:[0.5,0.6]}

#I need to read the freq table and store rsId as a key and frequency as value
all_sites={}
with open('%s' %(freq_tab) ,'r') as freq_file:
	next(freq_file)
	for c_row in freq_file:
		site=c_row.rstrip().split("\t")
		all_sites[site[1]]=float(site[4])

#We will read each frequency and assign a bin foreach rs
binned_sites=collections.defaultdict(list)
bin_size={}
exclude_sites={}

for bins in maf_bins.keys():
	for rs_id in all_sites.keys():
		if (all_sites[rs_id] >= maf_bins[bins][0] and all_sites[rs_id] < maf_bins[bins][1]):
			binned_sites[bins].append(rs_id)

#count how many sites we have in each bin and select 10% from each bin, picking random values
for bin_f in binned_sites.keys():
	bin_size[bin_f]=[len(binned_sites[bin_f]), int((len(binned_sites[bin_f])*20)/100)]
	#now we select randomly the 10% of sites to remove from each bin
	c_bin_sites=binned_sites[bin_f]
	random.shuffle(c_bin_sites)
	exclude_sites[bin_f]=c_bin_sites[:bin_size[bin_f][1]]

exclude_list=open('%s.exclude' %(freq_tab),'a')
for ex_bin in exclude_sites.keys():
	for ex_site in exclude_sites[ex_bin]:
		print >> exclude_list,'%s' %(ex_site)
