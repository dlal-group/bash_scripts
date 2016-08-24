#!/usr/bin/env python2.7
#script used to check overlap between two different legend files
import gzip 
import re 
import sys 
import subprocess as sub
import time
import collections
import itertools

#this script is meant to be usede to check stuff on the reference panel input files, like duplicates or triplicates row
# alleles mismatches or other stuff

# vcfdata = sys.stdin.readlines()
# vcfdata ="test_last_chr4_REF_keep.vcf"
start_time = time.time()
start_time1 = time.ctime(int(start_time))
print start_time1

legend1=sys.argv[1]
panel1=sys.argv[2]

legend2=sys.argv[3]
panel2=sys.argv[4]

chrom=sys.argv[5]
# chrom=1
outdir=sys.argv[6]
# outdir="/netapp02/data/imputation/INGI_TGP3/impute/FVG/MERGED/CLEANED"

# we need to read one legend file and transform all in a dictionary
legend1_dict={}

with gzip.open('%s' %(legend1) ,'r') as legend1_file:
	next(legend1_file)
	for c_row in legend1_file:
		site=c_row.rstrip().split(" ")
		legend1_dict[site.join("_")] = panel1

legend2_dict={}
with gzip.open('%s' %(legend2) ,'r') as legend2_file:
	next(legend2_file)
	for c_row2 in legend2_file:
		site2=c_row2.rstrip().split(" ")
		legend2_dict[site2.join("_")] = panel2

#now we define a dictionary in wich we will merge the two dict previously created
#so that we have the same key, but the value will be the name of one panel o the othe or both of them
all_panels=collections.defaultdict(list)

for k,v in chain(legend1.iteritems(), legend2.iteritems()):
	all_panels.setdefault(k, []).extend(v)

#now we count how many sites are in panel 1, how manin in panel 2 and how many in both, and print it
share_count = [(k, len(list(v))) for k, v in itertools.groupby(sorted(all_panels.values()))]

panel_share=open('%s/chr%s_%s_%s.share' %(outdir, chrom,panel1,panel2), 'w')
print >> panel_share,'%s %s %s ' %(panel_share[0][0],panel_share[1][0].panel_share[2][0])
print >> panel_share,'%s %s %s ' %(panel_share[0][1],panel_share[1][1].panel_share[2][1])
panel_share.close()

#now lets print the list of sites to keep			
print "Keeplist created!"
#now we need to read the stream from the vcf file and for each line
# decide if we want to keep it or not, based on matching fields

endtime1=time.time() - start_time

print time.ctime(int(endtime1))
