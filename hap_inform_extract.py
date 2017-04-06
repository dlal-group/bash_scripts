#!/usr/bin/env python
import gzip 
import re 
import sys 
import time
import collections
import itertools

input_file=sys.argv[1]
n_snp=sys.argv[2]
#input_file="ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.0.5.vcf.gz"
#n_snp=6


# lines= [if not re.match('#',i): i.split("\t") for i in current_file.readlines()]
current_file=gzip.open('%s' %(input_file), 'r')
lines=[]
for line in current_file:
	if not re.match('#',line):
		site=line.rstrip().split("\t")
		geno=site[9:]
		# need to take care of multiallelic sites
		geno=[w.replace('0', site[3]) for w in geno]
		if len(site[4]) > 1:
			multi=",".split(site[4])
			geno=[w.replace('1', multi[0]) for w in geno]
			geno=[w.replace('2', multi[1]) for w in geno]
		else
			geno=[w.replace('1', site[4]) for w in geno]
		# variant=[site[0],site[1],site[2], site[3],site[4],site[9:]]
		variant=[site[0],site[1],site[2], site[3],site[4],geno]
		lines.append(variant)

all_rs_list=[]
rs_list=[]
for i in range(0, len(lines)):
	j =i+int(n_snp)
	#rs_list=[]
	if j < len(lines):
		if int(lines[j][1]) -int(lines[i][1]) > 100: 
			if len(rs_list) != 0:
				# all_rs_list.append(list(set(rs_list)))
				rs_list.sort()
				all_rs_list.append(list(rs_list for rs_list,_ in itertools.groupby(rs_list)))
			rs_list=[]
			continue
		else: 
			for k in range(0, int(n_snp)):
				#print(lines[i+k][2], lines[i+k][4], lines[i+k][9])
				# rs_list.append(lines[i+k][2])
				rs_list.append(lines[i+k])
			while int(lines[j][1]) -int(lines[i][1]) <= 100:
				#print(lines[j][2], lines[j][4], lines[j][9])
				rs_list.append(lines[j])
				j +=1
				# rs_list.append(lines[j][2])
 
for snp_set in all_rs_list:
	print 'block\t%s' %(all_rs_list.index(snp_set))
	for snp in snp_set:
		print '%s\t%s\t%s\t%s\t%s\t%s' %(snp[0],snp[1],snp[2],snp[3],snp[4],"\t".join(snp[5]))

# print(all_rs_list)
