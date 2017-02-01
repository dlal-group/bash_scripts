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
import os

all_files=sys.argv[1:2]
chrom=sys.argv[3]
# chrom=10
# file3=sys.argv[3]
pop=sys.argv[4]
# pop="FVG"

#read the two info files, if the extension in gz, we need to read it as a gzipped file
# all_files=["/netapp02/data/imputation/INGI_TGP3/CARL/carl/MERGED/ALL/chr10.gen_info","/home/cocca/imputation/MERGED_INGI_TGP3_23012017/FVG/MERGED/ALL/chr10.gen_info"]
# info_file="/netapp02/data/imputation/INGI_TGP3/CARL/carl/MERGED/ALL/chr10.gen_info"
# to_merge_files={}
to_merge_files=collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
all_sites_keys=[]

sys.stderr.write(str(all_files))
sys.stderr.write(str(sys.argv))

for info_file in all_files:
	# info_file=all_files[0]
	pathfilename, file_extension = os.path.splitext(info_file)
	filename = os.path.basename(info_file)
	#we need to modify the filename since they can be the same for different imputation set
	key_filename="_".join([filename,str(all_files.index(info_file))])
	sys.stderr.write("Reading file " + filename + "...\n")
	if file_extension == ".gz" :
		current_file=gzip.open('%s' %(info_file), 'r')
		for line in current_file:
			if not re.match('snp_id',line):
				site=line.rstrip().split(" ")
				site_k=site[1]
				to_merge_files[site_k][key_filename] = {"chrom":chrom,"pos":site[2],"exp_a1_freq":site[5],"info":site[6],"type":site[8]}
	else :
		with open('%s' %(info_file) ,'r') as current_file:
			next(current_file)
			for line in current_file:
				site=line.rstrip().split(" ")
				site_k=site[1]
				to_merge_files[site_k][key_filename] = {"chrom":chrom,"pos":site[2],"exp_a1_freq":site[5],"info":site[6],"type":site[8]}

#create header:I need to adjust for the number of files selected for comparison
print 'CHROM\tPOS',
for info_file_name in all_files:
	filename = os.path.basename(info_file_name)
	#we need to modify the filename since they can be the same for different imputation set
	header_key_filename="_".join([filename,str(all_files.index(info_file_name))])
	print '\t%s\t%s\t%s' %("info_"+header_key_filename, "type_"+header_key_filename, "exp_a1_af_"+header_key_filename),

for intersect_site in to_merge_files.keys():
	if len(to_merge_files[intersect_site]) > 1:
		f_k=to_merge_files[intersect_site].keys()[0]
		chrom=to_merge_files[intersect_site][f_k]['chrom']
		pos=to_merge_files[intersect_site][f_k]['pos']
		print '%s\t%s' %(chrom,pos),
		for info_filename in to_merge_files[intersect_site].keys():
			v_info=to_merge_files[intersect_site][info_filename]['info']
			v_type=to_merge_files[intersect_site][info_filename]['type']
			a1_freq=to_merge_files[intersect_site][info_filename]['exp_a1_freq']
			print '\t%s\t%s\t%s' %(v_info,v_type,a1_freq),
		print '\n'
