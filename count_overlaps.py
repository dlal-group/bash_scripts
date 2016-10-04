#!/usr/bin/env python2.7
#Script used to summarize data from bcftools isec
import gzip 
import re 
import sys 
import subprocess as sub
import time
import collections
import itertools


#we want to count overlaps from the esults written in files from bcftools isec
# need to manage different population numbers

start_time = time.time()
start_time1 = time.ctime(int(start_time))
print start_time1

isec_sites_file=sys.argv[1]
# isec_sites_file="/lustre/scratch113/projects/esgi-vbseq/09022016_PANEL_SESOURCES/INGI/UNION/1/sites_test.txt"
out_file=sys.argv[2]
# out_file="/lustre/scratch113/projects/esgi-vbseq/09022016_PANEL_SESOURCES/INGI/UNION/1/count.txt"
n_pop=sys.argv[3]
# n_pop=3

#we need to use this to get out the correct header
pops=sys.argv[4:-1]
# pops=["CARL","VBI","FVG"]


# first, we need to create a dictionary with position as the key, than splitted by 
# variant type: well add all the data in the value

all_sites=collections.defaultdict(lambda: collections.defaultdict(list))
with open('%s' %(isec_sites_file) ,'r') as isec_file:
	for c_row in isec_file:
		site=c_row.rstrip().split("\t")
		if len(site[2]) != len(site[3]):
			all_sites[site[1]]["INDEL"].append([site[0],site[1],site[2],site[3],site[4]])
		else:
			all_sites[site[1]]["SNP"].append([site[0],site[1],site[2],site[3],site[4]])

#1)Based on the number of populations we need to count how many items we have 
# shared or private among the populations

#we need to define the population code based on population number
#need a dictionary!
pops_share={}
pops_comb=[]
for L in range(1, len(pops)+1):
	for subset in itertools.combinations(pops, L):
		pops_comb.append(list(subset))

for comb in pops_comb:
	comb_share=["0"] * n_pop
	if len(comb) == 1:
		pop_i=pops.index(comb[0])
		comb_share[pop_i] = "1"
		pops_share[pops[pop_i]]="".join(comb_share)
	else:
		for pop in comb:
			pop_i=pops.index(pop)
			comb_share[pop_i] = "1"
		pops_share["_".join(comb)]="".join(comb_share)	

#invert keys and values for a better check
pops_share_inv = {y:x for x,y in pops_share.iteritems()}

#now we need to count what we have and assign everything to the correct category
all_count={}
all_count["SNP"]={y:0 for x,y in pops_share.iteritems()}
all_count["INDEL"]={y:0 for x,y in pops_share.iteritems()}

for key in all_sites:
	for v_type in all_sites[key]:
		dupe_map=[]
		dupe_sum=[]
		if len(all_sites[key][v_type]) == 1:
			all_count[v_type][all_sites[key][v_type][0][4]] = all_count[v_type][all_sites[key][v_type][0][4]]+1 
		else:
			#we need to manage the duplicates
			for dupe in all_sites[key][v_type]:
				dupe_map.append(list(bool(int(x)) for x in list(dupe[4])))
			for dupe_i in range(0,n_pop):
				dupe_b_sum=[]
				for dupe_b in dupe_map:
					dupe_b_sum.append(dupe_b[dupe_i])
				dupe_sum.append(dupe_b_sum)
			dupe_pop_count="".join(list(str(int(bool(sum(x)))) for x in dupe_sum))
			all_count[v_type][dupe_pop_count] = all_count[v_type][dupe_pop_count]+1 

#now we need to write the outfile
for vtype in all_count:
	share_count=open('%s.%s.tab' %(out_file,vtype), 'w')
	#match count key with pop
	for p_key in all_count[vtype]:
		print >> share_count,'%s %s' %(pops_share_inv[p_key],all_count[vtype][p_key])
	share_count.close()


#now lets print the list of sites to keep			
print "Keeplist created!"
#now we need to read the stream from the vcf file and for each line
# decide if we want to keep it or not, based on matching fields

endtime1=time.time() - start_time

print time.ctime(int(endtime1))
