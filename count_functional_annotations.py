import sys 
import gzip
import re 
import os 

chr=sys.argv[1]
#tags=['REGULATORY_REGION','UPSTREAM','WITHIN_NON_CODING_GENE','DOWNSTREAM','INTRONIC','INTERGENIC','NON_SYNONYMOUS_CODING','SYNONYMOUS_CODING','5PRIME_UTR','3PRIME_UTR','NMD_TRANSCRIPT','SPLICE_SITE','STOP_GAINED','ESSENTIAL_SPLICE_SITE','WITHIN_MATURE_miRNA','STOP_LOST','PARTIAL_CODON','CODING_UNKNOWN']
tags=['WITHIN_NON_CODING_GENE','INTRONIC','DOWNSTREAM','UPSTREAM','INTERGENIC','REGULATORY_REGION','NON_SYNONYMOUS_CODING','5PRIME_UTR','SYNONYMOUS_CODING','3PRIME_UTR','NMD_TRANSCRIPT','ESSENTIAL_SPLICE_SITE','SPLICE_SITE','STOP_GAINED','STOP_LOST','WITHIN_MATURE_miRNA']
diktag_snp={}
diktag_indel={}
for cat in tags: 
	diktag_snp[cat]=0
	diktag_indel[cat]=0

#this work by chr!!
print chr 
#pat tho vcf 
inputfile=gzip.open('/nfs/users/nfs_m/mc14/lustre110_home/GENOTIPI/COMPARISON/NOT_OVERLAPPING/PUTATIVE_NOVEL/NEW_RUN/former_rsID_filtered/UK10K_FILTERED/gte2_indivs/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr%s.re_ann.NOT_OVERLAP.NO_RSID.no_UK10K.gte2_indivs.vcf.gz' %(chr), 'r')

for line in inputfile:
	x=line.split() 
	if x[0]==chr:
		#print line  
		for cat in tags:
			if re.search(cat, line): 
				if re.search('INDEL', line): diktag_indel[cat]+=1
				else: diktag_snp[cat]+=1
out=open('/nfs/users/nfs_m/mc14/lustre110_home/GENOTIPI/COMPARISON/NOT_OVERLAPPING/PUTATIVE_NOVEL/NEW_RUN/former_rsID_filtered/UK10K_FILTERED/gte2_indivs/FUNKANN/funkann.chr%s.table' %(chr), 'w')
sys.stdout=out
print 'CAT\tINDEL\tSNP' 
					
for cat in tags: print '%s\t%s\t%s' %(cat,diktag_indel[cat], diktag_snp[cat])

