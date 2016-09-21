import sys
import gzip
import re
import os

chr=sys.argv[1]
#first_path=sys.argv[2]
#second_path=sys.argv[3]


sitivariabilipop=[]
sitivariabilivbi=[]

if re.match('23',chr):
	chr='X'

my_POPfile=open('/lustre/scratch109/sanger/mc14/GENOTIPI/COMPARISON/VBSEQ_QC/CEU/CEU.chr%s.variable.site' %(chr), 'r')

my_VBIfile=gzip.open('/lustre/scratch110/sanger/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr%s.re_ann.vcf.gz' %(chr), 'r')

for line in my_POPfile:
	x=line.split()
	if not re.match('#',line):
		sitivariabilipop.append(x[1])


for line in my_VBIfile:
	x=line.split()
	if not re.match('#',line):
		sitivariabilivbi.append(x[1])


matches=set(sitivariabilipop) & set(sitivariabilivbi)

pop_private=set(sitivariabilipop) - set(sitivariabilivbi)


private=open('/lustre/scratch110/sanger/mc14/GENOTIPI/COMPARISON/NOT_OVERLAPPING/VBI_vs_CEU/TAB/private_ceu_sites.chr%s.list' %(chr), 'w')
sys.stdout=private
print 'CHR\tPOS'

for pos in pop_private: print '%s\t%s' %(chr,pos)


overlap=open('/lustre/scratch110/sanger/mc14/GENOTIPI/COMPARISON/OVERLAPPING/VBI_vs_CEU/TAB/overlapping_sites.chr%s.list' %(chr), 'w')
sys.stdout=overlap
print 'CHR\tPOS'

for pos in matches: print '%s\t%s' %(chr,pos)
