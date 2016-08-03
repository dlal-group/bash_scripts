#!/software/bin/python2.7
#script used to re filter vcf files for reference generation
import gzip 
import re 
import sys 
import subprocess as sub
import time


#this script is meant to be usede to check stuff on the reference panel input files, like duplicates or triplicates row
# alleles mismatches or other stuff

# vcfdata = sys.stdin.readlines()
# vcfdata ="test_last_chr4_REF_keep.vcf"
file_name=sys.argv[1]
chrom=sys.argv[2]
pop=sys.argv[3]
mode=sys.argv[4]

#we need to work in different ways depending on the MODE parameter
#it could have two values => INFO or GEN
if mode == "INFO":
	
	
elif mode == "GEN":




#first thing we want to check for duplicates/triplicates rows
# the input data has to be in the form of a table with:
# CHROM POS REF ALT INFO/DP4 INFO/DP4 INFO/DP INFO/HOB INFO/ICB INFO/IDV
all_list_name=all_list.split("/")[-1]

all_sites_dict={}
# pos_dict={}

start_time1 = time.time()

print start_time1

with open('%s' %(all_list) , 'r') as keep_file:
	for keep_row in keep_file:
		keep_line=keep_row.rstrip().split("\t")yep, mee
		# check_string=single_line[0]+single_line[1]+single_line[4]+single_line[5]+single_line[6]+single_line[7]+single_line[8]+single_line[9]
		sites_key=(keep_line[0],keep_line[1],keep_line[2],keep_line[3],keep_line[4])
		# first mark all duplicates by position
		all_sites_dict[sites_key] = keep_line[5]

endtime1=time.time() - start_time1
print endtime1

print "Dictionaries created!"
#now we need to read the stream from the vcf file and for each line
# decide if we want to keep it or not, based on matching fields
keep_in=open('%s/%s.%s.%s.to_keep.vcf' %(outdir, all_list_name, cohort, mode), 'w')
# keep_out=open('%s/%s.%s.%s.to_remove.vcf' %(outdir, all_list_name, cohort, mode), 'w')
# keep_in=open('%s.keep_TEST.vcf' %(all_list_name), 'w')

# with open('%s' %(vcfdata) , 'r') as vcfdata_file:
	# for vcf_row in vcfdata_file:
for vcf_row in vcfdata:
		# infile = open('%s' %(vcfdata) , 'r')
		# vcf_row = infile.readline()
	vcf_line=vcf_row.rstrip().split("\t")
	info_field=vcf_line[7].rstrip().split(";")
	# info_key=";".join([info_field[0],info_field[1],info_field[2],info_field[3],info_field[4]])
	info_key=";".join([info_field[0],info_field[1]])
	sites_key_vcf=(vcf_line[0],vcf_line[1],vcf_line[3],vcf_line[4],info_key)
	if sites_key_vcf in all_sites_dict:
		print >> keep_in,'%s' %(vcf_row.rstrip())
	# else :
	# 	print >> keep_out,'%s' %(vcf_row.rstrip())


keep_in.close()