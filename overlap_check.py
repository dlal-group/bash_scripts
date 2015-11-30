#!/software/bin/python2.7
#script used to extract variants overlapping in different populations after bcftools isec command
import gzip 
import re 
import sys 
import subprocess as sub


#this version works on a single chromosome, so we can submit a job array
# var_list="/lustre/scratch113/projects/esgi-vbseq/27112015_INGI_REF_PANEL/VBI/22.vcf.gz.snps_ac1dp5.tab"
# overlap_list="/lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/INGI/UNION/22/sites.txt"
# outdir="/lustre/scratch113/projects/esgi-vbseq/27112015_INGI_REF_PANEL/VBI"

cohort=sys.argv[1]
var_list=sys.argv[2]
overlap_list=sys.argv[3]
outdir=sys.argv[4]
var_list_name=var_list.split("/")[-1]

#read overlap list and create a dictionary with chr, pos and score, where score is the sum of the elements in the
# last column
over_dict={}
for row in open('%s' %(overlap_list) , 'r'):
	over_line=row.rstrip().split("\t")
	over_dict[(over_line[0],over_line[1])] = sum(int(x) for x in over_line[4] if x.isdigit())

#now we need to write 2 files:
#1) a tabbed file with overlapping variants
#2) a tabbed file with non overlapping variants wich we will check against other stuff
over_out=open('%s/%s.%s.over.tab' %(outdir, var_list_name, cohort), 'w')
not_over_out=open('%s/%s.%s.not_over.tab' %(outdir, var_list_name, cohort), 'w')


over_var={}
not_over_var=[]
for row in open('%s' %(var_list) , 'r'):
	var_line=row.rstrip().split("\t")
	if (over_dict[(var_line[0],var_line[1])] > 1):
		over_var[(var_line[0],var_line[1])]=over_dict[(var_line[0],var_line[1])]
		print >> over_out,'%s\t%s' %(var_line[0],var_line[1])
	else:
		not_over_var.append((var_line[0],var_line[1]))
		print >> not_over_out,'%s\t%s' %(var_line[0],var_line[1])
