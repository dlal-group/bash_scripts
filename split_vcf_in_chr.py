import sys
import re
import gzip

chr=sys.argv[1]
out=gzip.open('/nfs/users/nfs_m/mc14/lustre_home/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/SPLITTED_VCF/chr%s.esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz' %( chr), 'w')
sys.stdout=out 

for line in gzip.open('/nfs/users/nfs_m/mc14/lustre_home/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.vcf.gz',  'r'):
	strline=line.rstrip()
	if re.search('#' , line ): print strline
	elif re.match('%s\t' %(chr), line ): print strline
