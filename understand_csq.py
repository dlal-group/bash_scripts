#!/software/bin/python2.7

import gzip 
import re 
import sys 

listaconseq=[]

for chr in range (1, 23): 
	print chr 
	# for line in gzip.open('/nfs/users/nfs_m/mc14/lustre110_home/GENOTIPI/COMPARISON/NOT_OVERLAPPING/PUTATIVE_NOVEL/NEW_RUN/former_rsID_filtered/UK10K_FILTERED/gte2_indivs/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr%s.re_ann.NOT_OVERLAP.NO_RSID.no_UK10K.gte2_indivs.vcf.gz' %(chr) , 'r'): 
	for line in gzip.open('/lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/FIXED_ALT/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr%s.re_ann.alt_fixed.vcf.gz' %(chr) , 'r'): 
		if re.match('\d+\t', line):
			x=line.split('\t')
			y=x[7].split(';')
			for item in y:
				if re.match('CSQ', item):
					#print item  
					z=item.split('+')
					for conseq in z:
						if not re.search('GERP', conseq):  
							w=conseq.split(':')
							if re.search(',', w[2]) : 
								a=w[2].split(',')
								for tipi in a: 
									if not tipi in listaconseq: listaconseq.append(tipi) 
							else: 
								if not w[2] in listaconseq: listaconseq.append(w[2])

out=open('/lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/FIXED_ALT/consequences.list', 'w')
sys.stdout=out
for csq in listaconseq: print csq
						 
				
