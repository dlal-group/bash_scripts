#!/usr/bin/env python2.7
#annotate variants after genabel analyses using 1000GPh3 data
import argparse
import gzip 
import re 
import sys 
import subprocess as sub
import time
import collections
import itertools
from itertools import chain
from datetime import timedelta
import select
from subprocess import call
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool 


# args=parser.parse_args()
roh_file=sys.argv[1] #roh_file="/lustre/scratch113/projects/carl_seq/08072016_paperONE/FILTERED_MAF/18112016_IBDSEQ_r20.15_mac1_CEU/CHR1/CEU.roh.hbd"
vcf_file=sys.argv[2] #vcf_file="/lustre/scratch113/projects/carl_seq/08072016_paperONE/TGP3_INGI_MERGE/INGI_TGP3_filt.vcf.gz"
filtered_snps=sys.argv[3] #filtered_snps="/lustre/scratch113/projects/carl_seq/08072016_paperONE/FILTERED_MAF/18112016_IBDSEQ_r20.15_mac1_CEU/CHR1/CEU.roh.r2.filtered.gz"
out_file=sys.argv[4] #out_file="/lustre/scratch113/projects/carl_seq/08072016_paperONE/FILTERED_MAF/18112016_IBDSEQ_r20.15_mac1_CEU/CHR1/CEU.roh.hbd.extended"

######################################## CLASS and methods definitions ###################################################
#define a class for ROH segments
class ROH(object):
	chrom=''
	start=''
	end=''
	lod=''
	hets=''
	filtered_hets=''
	density=''
	filtered_density=''
	"""docstring for ROH"""
	def __init__(self,chrom,start,end,lod):
		super(ROH, self).__init__()
		self.chrom = chrom
		self.start = start
		self.end = end
		self.lod = lod
		self.length= int(end) - int(start) + 1
		self.region= "".join([str(chrom),":",str(start),"-",str(end)])
		self.hets= None
		self.filtered_hets= None
		self.density= None
		self.filtered_density= None

def make_roh(chrom,start,end,lod):
    roh = ROH(chrom,start,end,lod)
    return roh

def ann_reader(d,ann_line):
	x=ann_line.split("\t")
	if not re.match('#',line):
		ann_key="_".join([x[0],x[1],x[3],x[4]])
		d[ann_key] = x[2]

def check_density(s,vcf,filter_list=None):
	#filter is a list of sites' positions to remove,if specified
	var_count_cmd=" ".join(["bcftools view -H -G -r",s.region,str(vcf),"| wc -l"])
	bcft=sub.Popen(var_count_cmd,shell=True,stdout=sub.PIPE)
	sites=bcft.communicate()
	if filter_list:
		#get the sites from the filter list that are in my region 
		#remove them from the count
		excluded_count_cmd="".join(["zcat ",str(filter_list)," | awk '$2>=",s.start," && $2<=",s.end,"'| wc -l"])
		exclude_count=sub.Popen(excluded_count_cmd,shell=True,stdout=sub.PIPE)
		exclude_sites=exclude_count.communicate()
		filtered_density = int(sites[0].strip()) - int(exclude_sites[0].strip())
		density = sites[0].strip()
	else:
		density = sites[0].strip()
		filtered_density = sites[0].strip()
	return density, filtered_density


def check_het(sample,s,vcf,filter_list=None):
	#filter is a list of sites' positions to remove,if specified
	gen_count_cmd=" ".join(['bcftools query -f"[%GT]\n" -s',sample,'-r',s.region,vcf])
	gen_bcft=sub.Popen(gen_count_cmd,shell=True,stdout=sub.PIPE)
	s_geno=gen_bcft.communicate()
	hets=gen_conv(s_geno[0].split()).count(1)
	if filter_list:
		#get the sites from the filter list that are in my region 
		#remove them from the count
		# excluded_count_cmd="".join(["zcat ",str(filter_list)," | awk '$2>=",s.start," && $2<=",s.end,"'| wc -l"])
		# exclude_count=sub.Popen(excluded_count_cmd,shell=True,stdout=sub.PIPE)
		# exclude_sites=exclude_count.communicate()
		filtered_hets = ''
	else:
		filtered_hets = hets
	return hets, filtered_hets

#4) Add also het number (before and after)
def gen_conv(genos):
	return [ sum(map(int,geno.split("|"))) for geno in genos]
	
#########################################################################################
#1) Read exclusion file
sys.stderr.write('Reading exclusion file...\n')
start_time = time.time()

# pool=ThreadPool(int(threads))
current_filtered=gzip.open('%s' %(filtered_snps), 'r')
all_filtered={}

for line in current_filtered:
	x=line.split()
	ann_key="_".join([x[0],x[1],x[3],x[4]])
	all_filtered[ann_key] = x[2]

elapsed_time = time.time() - start_time
sys.stderr.write('exclusion file read in '+ str(timedelta(seconds=elapsed_time)) +'...\n')
##########################################################################################
#2) Read roh file
sys.stderr.write('Reading ROH file...\n')
start_time = time.time()

all_ROH=collections.defaultdict(lambda: collections.defaultdict(dict))
with open('%s' %(roh_file) ,'r') as current_roh:
	for line in current_roh:
		x=line.split()
		sample=x[0]
		roh_key="_".join([x[4],x[5],x[6]])
		all_ROH[sample][roh_key] = make_roh(x[4],x[5],x[6],x[7])


elapsed_time = time.time() - start_time
sys.stderr.write('ROH file read in '+ str(timedelta(seconds=elapsed_time)) +'...\n')

#########################################################################################
#3) Count number of sites before and after removing the filtered site
#4) Print the in the output file

outfile=open(out_file,'w')

sys.stderr.write('Create extended ROH file...\n')
start_time = time.time()

for sample in all_ROH:
	for roh in all_ROH[sample]:
		all_ROH[sample][roh].density,all_ROH[sample][roh].filtered_density = check_density(all_ROH[sample][roh],vcf_file,filtered_snps)
		all_ROH[sample][roh].hets,all_ROH[sample][roh].filtered_hets=check_het(sample,all_ROH[sample][roh],vcf_file,filtered_snps)
		print >> outfile, '%s\t1\t%s\t2\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(sample,sample,all_ROH[sample][roh].chrom,all_ROH[sample][roh].start,all_ROH[sample][roh].end,all_ROH[sample][roh].lod,all_ROH[sample][roh].density,all_ROH[sample][roh].filtered_density,all_ROH[sample][roh].hets)

elapsed_time = time.time() - start_time
sys.stderr.write('Extended ROH file created in '+ str(timedelta(seconds=elapsed_time)) +'...\n')
