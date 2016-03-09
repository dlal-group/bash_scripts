#!/software/bin/python2.7
#script used to extract variants overlapping in different populations after bcftools isec command
import gzip 
import re 
import sys 
import subprocess as sub


#this script is meant to be usede to check stuff on the reference panel input files, like duplicates or triplicates row
# alleles mismatches or other stuff

# data = sys.stdin.read()
cohort=sys.argv[1]
all_list=sys.argv[2]
dupe_pos_list=sys.argv[3]
dupe_sites_list=sys.argv[4]
# all_list="test.tab"
outdir=sys.argv[5]
mode=sys.argv[6]

#first thing we want to check for duplicates/triplicates rows
# the input data has to be in the form of a table with:
# CHROM POS REF ALT INFO/DP4 INFO/DP4 INFO/DP INFO/HOB INFO/ICB INFO/IDV
all_list_name=all_list.split("/")[-1]
sites_dict={}
pos_dict={}

for pos_row in open('%s' %(dupe_pos_list) , 'r'):
	pos_line=pos_row.rstrip().split("\t")
	# check_string=single_line[0]+single_line[1]+single_line[4]+single_line[5]+single_line[6]+single_line[7]+single_line[8]+single_line[9]
	# first mark all duplicates by position
	pos_dict[pos_line[0]] = pos_line[1]

for dupe_row in open('%s' %(dupe_sites_list) , 'r'):
	dupe_line=dupe_row.rstrip().split("\t")
	# check_string=single_line[0]+single_line[1]+single_line[4]+single_line[5]+single_line[6]+single_line[7]+single_line[8]+single_line[9]
	# first mark all duplicates by position
	site_key=(dupe_line[0],dupe_line[1],dupe_line[4],dupe_line[5],dupe_line[6],dupe_line[7],dupe_line[8])
	sites_dict[site_key] = dupe_line[9]

print "Dictionaries created!"
#we need to write 2 files:
#1) a tabbed file with keep variants
#2) a tabbed file with excluded variants
keep_out=open('%s/%s.%s.%s.to_keep.tab' %(outdir, all_list_name, cohort, mode), 'w')
# keep_out=open('/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/VBI/PANEL/%s.keep.tab' %(all_list_name), 'w')
#all those with 0 positions duplicates, we want to keep them
#plus, we want to keep all those duplicates by positions which also have a duplicate by row
for row in open('%s' %(all_list) , 'r'):
	# infile = open('%s' %(all_list) , 'r')
	# row = infile.readline()
	var_line=row.rstrip().split("\t")
	pos_key=var_line[1]
	# for pos_key in pos_dict.keys():
		# pos_key=pos_dict.keys()[0]
	contained = [key for key in sites_dict.keys() if pos_key in key]
	# for sites_key in sites_dict.keys():
	for sites_key in contained:
		# sites_key=contained[0]
		if pos_key in sites_key:
			# we need to create 2 sets for the intersection
			s1=set(sites_key)
			s2=set(var_line)
			# if this intersection is equal to the line I'm reading and to the site key, than go on
			if s1.intersection(s2) == s1:
			# if sites_key in var_line:
				if pos_dict[pos_key] == 0:
				#need to keep the site since there are not position duplicates
				# ['1', '231955500', 'A', 'AT', '1162,1312,57,90', '2839', '0.00578704', '0.906249', '9']
					print >> keep_out,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(var_line[0],var_line[1],var_line[2],var_line[3],var_line[4],var_line[5],var_line[6],var_line[7],var_line[8],"KEEP")
				else :
					#if we have dupolicates by position , we need to check them
					if sites_dict[sites_key] != 0:
						# we'll keep the duplicates only if they're from the same site (splitted multiallelic site)
						print >> keep_out,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(var_line[0],var_line[1],var_line[2],var_line[3],var_line[4],var_line[5],var_line[6],var_line[7],var_line[8],"KEEP")

keep_out.close()
