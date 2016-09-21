#!/software/bin/python2.7
#script to define chr_chunks from a vcf file of input
import io
import gzip 
import re 
import sys 
import subprocess as sub
from subprocess import Popen, PIPE 


#it works on a chr base or on a all genome vcf

# cohort="VBI"
# var_list="/lustre/scratch113/projects/esgi-vbseq/27112015_INGI_REF_PANEL/VBI/22.vcf.gz.snp_ac1dp5.tab"
# overlap_list="/lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/INGI/UNION/22/sites.txt"
# outdir="/lustre/scratch113/projects/esgi-vbseq/27112015_INGI_REF_PANEL/VBI"
# mode="snp"
# mode="indel"
vcf=sys.argv[1]

chunk_size = 5000000 #<- We wan chunks of 5 Mb ...

for i_line in gzip.open('%s' %(vcf) , 'r'): 
	if not re.match('^#', i_line):
		# I just want the first position
		x=i_line.split('\t')
		chr=x[0]
		i_pos=int(x[1])
		break

# Implement tabix <file> <chr> | tail -1 
# head = sub.Popen(['head', '-1'], stdin=gzip.stdout, stdout=PIPE) 
gzip = sub.Popen(['tabix', vcf, chr], stdout=PIPE) 
tail = sub.Popen(['tail', '-1'], stdin=gzip.stdout, stdout=PIPE) 
f_line = tail.communicate()[0] 
f_x=f_line.split('\t')
f_pos=int(f_x[1])

chunks_file=open('%s.chunks' %(vcf), 'w')

#now write chunks accordingly
for region in range(i_pos,f_pos,chunk_size):
	if (region+chunk_size <= f_pos):
		print >> chunks_file,'%s:%s-%s' %(chr,region,region+chunk_size-1)
	else:
		print >> chunks_file,'%s:%s-%s' %(chr,region,f_pos)
