#!/software/bin/python2.7
#script used to extract consequence annotation from vcf files. From Enza Colonna and moddified by Max Cocca.
import gzip 
import re 
import sys 
import subprocess as sub

#this version works on a single chromosome, so we can submit a job array
in_path=sys.argv[1]
# in_path="/lustre/scratch113/projects/esgi-vbseq/08092015/variant_calling/lists/chr20-pooled.list"

for line in open('%s' %(in_path) , 'r'):
	bam_file=line.rstrip()
	# cmd=['/software/hgi/pkglocal/samtools-1.2/bin/samtools','view','-H',bam_file,'|','grep','^@RG','|','head','-1']
	cmd=['/software/hgi/pkglocal/samtools-1.2/bin/samtools','view','-H',bam_file]
	p = sub.Popen(cmd,stdout=sub.PIPE,stderr=sub.PIPE)
	output, errors = p.communicate()
	for tag in output.rstrip().split("\n"):
		c_tag=tag.rstrip()
		if re.match('@RG', c_tag):
			x=c_tag.split('\t')
			for bam_h in x:
				if re.match('SM:', bam_h):
					clinic=bam_h.split(":")
					sample=clinic[1]
					print bam_file + " " + sample
			break
			