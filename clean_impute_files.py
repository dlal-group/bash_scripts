#!/usr/bin/env python2.7
#script used to re filter vcf files for reference generation
import gzip 
import re 
import sys 
import subprocess as sub
import time
import collections


#this script is meant to be usede to check stuff on the reference panel input files, like duplicates or triplicates row
# alleles mismatches or other stuff

# vcfdata = sys.stdin.readlines()
# vcfdata ="test_last_chr4_REF_keep.vcf"
start_time = time.time()
start_time1 = time.ctime(int(start_time))
print start_time1

file_prefix=sys.argv[1]
# file_prefix="chr1.62.gen"
chrom=sys.argv[2]
# chrom=1

pop=sys.argv[3]
# pop="FVG"
outdir=sys.argv[4]
# outdir="/netapp02/data/imputation/INGI_TGP3/impute/FVG/MERGED/CLEANED"

# first, we need to create a dictionary with chr_pos_a1_a2 as key and [chr,pos,a1,a2,variant type] as value
# this is the same for both info and gen files, but we'll have a gzipped file
file_info=file_prefix + "_info"
file_gen=file_prefix + ".gz"

all_sites=collections.defaultdict(lambda: collections.defaultdict(list))
with open('%s' %(file_info) ,'r') as info_file:
	next(info_file)
	for c_row in info_file:
		site=c_row.rstrip().split(" ")
		if len(site[3]) != len(site[4]):
			all_sites[site[2]]["INDEL"].append([chrom,site[2],site[3],site[4]])
		else:
			all_sites[site[2]]["SNP"].append([chrom,site[2],site[3],site[4]])

	#1)We need to remove duplicates for INDELs and SNPs
	# we'll print al list of sites TO KEEP formatted as:
	# CHR POS REF ALT
	# first we need to check all position with both snp and indel
keep_in=open('%s/chr%s.gen.to_keep' %(outdir, chrom), 'w')
keep_variants = []
for key in all_sites:
	if len(all_sites[key]) > 1:
		#this means we have both snps and indels in our dictionary's entry
		if len(all_sites[key]['SNP']) == 1:
			variant = map(str,all_sites[key]['SNP'][0])
			# keep_variants.append(variant)
			print >> keep_in,'%s' %(' '.join(variant))
			# print variant[0],' ',variant[1],' ',variant[2],' ',variant[3],' ','SNP'
	else:
		#we only have a single type of variant in that position
		for v_type in all_sites[key]:
			#we keep it if its not multiallelic
			if len(all_sites[key][v_type]) == 1:
				variant = map(str,all_sites[key][v_type][0])
				# keep_variants.append(variant)
				print >> keep_in,'%s' %(' '.join(variant))
				# print variant[0],' ',variant[1],' ',variant[2],' ',variant[3],' ',v_type

#now lets print the list of sites to keep			
print "Keeplist created!"
#now we need to read the stream from the vcf file and for each line
# decide if we want to keep it or not, based on matching fields

endtime1=time.time() - start_time

print time.ctime(int(endtime1))

keep_in.close()