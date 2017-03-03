#!/usr/bin/env python
import gzip 
import re 
import sys 
import time
import collections

input_file=sys.argv[1]
# input_file="ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.0.5.vcf.gz"


current_file=gzip.open('%s' %(input_file), 'r')
# lines= [if not re.match('#',i): i.split("\t") for i in current_file.readlines()]
lines=[]
for line in current_file:
	if not re.match('#',line):
		site=line.rstrip().split("\t")
		variant=[site[0],site[1],site[2]]
		lines.append(variant)

all_rs_list=[]
rs_list=[]
for i in range(0, len(lines)):
	j =i+6
	#rs_list=[]
	if j < len(lines):
		if int(lines[j][1]) -int(lines[i][1]) > 100: 
			if len(rs_list) != 0:
				all_rs_list.append(list(set(rs_list)))
			rs_list=[]
			continue
		else: 
			for k in range(0, 6):
				#print(lines[i+k][2], lines[i+k][4], lines[i+k][9])
				# rs_list.append(lines[i+k][2])
				rs_list.append(lines[i+k])
			while int(lines[j][1]) -int(lines[i][1]) <= 100:
				#print(lines[j][2], lines[j][4], lines[j][9])
				j +=1
				# rs_list.append(lines[j][2])
				rs_list.append(lines[j])
 
for snp_set in all_rs_list:
	print '%s\t%s\t%s\t%s' %(snp_set[0],),

print(all_rs_list)