#!/usr/bin/env python2.7
#script used to re filter vcf files for reference generation
import gzip 
import re 
import sys 
import time
import collections
import itertools
# import enumerate


#prepare data in tab format to be used to add annotations in VCF files
#right now we have six colunms TODO: extend to different number of cols

anno_file = sys.stdin.readlines()
ref_alleles=sys.argv[1]
# ref_alleles="1.MULTI_SNP.list"
# annotdata ="1.UNFORMATTED.CADD.tab"

# Easy peasy:create a dictionary with position and ref alelle as key
# Values will be annotations

# THIS IS MEANT to read the multiallelic table used to sort the alleles
all_multi_ref=collections.defaultdict((list))
with open('%s' %(ref_alleles) ,'r') as ref_file:
	for r_row in ref_file:
		r_site=r_row.rstrip().split("\t")
		all_multi_ref[r_site[1]] = [r_site[0],r_site[1],r_site[2],r_site[3]]

all_sites=collections.defaultdict(lambda: collections.defaultdict(list))
# FOR TESTING POURPUSES ONLY...or to use file name as argument
# with open('%s' %(annotdata) ,'r') as anno_file:
# 	next(anno_file)
# 	for c_row in anno_file:
# 		site=c_row.rstrip().split(" ")
# 		if len(site[2]) != len(site[3]):
# 			if all_multi_ref[site[1]]:
# 				all_sites[site[1]]["INDEL"].append([site[0],site[1],site[2],site[3],site[4],site[5],all_multi_ref[site[1]][3].split(",").index(site[3])])
# 			else:
# 				all_sites[site[1]]["INDEL"].append([site[0],site[1],site[2],site[3],site[4],site[5]])
# 		else:
# 			if all_multi_ref[site[1]]:
# 				all_sites[site[1]]["SNP"].append([site[0],site[1],site[2],site[3],site[4],site[5],all_multi_ref[site[1]][3].split(",").index(site[3])])
# 			else:
# 				all_sites[site[1]]["SNP"].append([site[0],site[1],site[2],site[3],site[4],site[5]])



#loop to leave in to use stdandard input 
for linenum, c_row in enumerate(anno_file):
	if linenum != 0:
		site=c_row.rstrip().split(" ")
		if len(site[2]) != len(site[3]):
			if all_multi_ref[site[1]]:
				all_sites[site[1]]["INDEL"].append([site[0],site[1],site[2],site[3],site[4],site[5],all_multi_ref[site[1]][3].split(",").index(site[3])])
			else:
				all_sites[site[1]]["INDEL"].append([site[0],site[1],site[2],site[3],site[4],site[5]])
		else:
			if all_multi_ref[site[1]]:
				all_sites[site[1]]["SNP"].append([site[0],site[1],site[2],site[3],site[4],site[5],all_multi_ref[site[1]][3].split(",").index(site[3])])
			else:
				all_sites[site[1]]["SNP"].append([site[0],site[1],site[2],site[3],site[4],site[5]])

#We need to collapse overlapping annotations
for key in all_sites:
# key='230571673'
	for v_type in all_sites[key]:
		# v_type='SNP'
		#we need to behave in a different way if it's multiallelic
		if len(all_sites[key][v_type]) == 1:
			variant = map(str,all_sites[key][v_type][0])
			print >> sys.stdout , '%s' %("\t".join(variant))
			# print variant[0],' ',variant[1],' ',variant[2],' ',variant[3],' ',v_type
		else:
			#we need to collapse stuff and remove duplicates and print, all sorted by the index value
			(all_sites[key][v_type]).sort(key=lambda x: x[6])
			collapsed = zip(*all_sites[key][v_type])
			collapsed_all=[]
			#we don't want to carry on the index used for the sorting
			for a_num, a in enumerate(collapsed[:-1]):
				if a_num < 3:
					collapsed_all.append(list(set(a)))
				else:
					collapsed_all.append(list(a))
			print >> sys.stdout, '%s' %("\t".join([",".join(sublist) for sublist in collapsed_all]))
