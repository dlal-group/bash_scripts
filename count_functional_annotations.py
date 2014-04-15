import sys 
import gzip
import re 
import os 

chr=sys.argv[1]
#tags=['REGULATORY_REGION','UPSTREAM','WITHIN_NON_CODING_GENE','DOWNSTREAM','INTRONIC','INTERGENIC','NON_SYNONYMOUS_CODING','SYNONYMOUS_CODING','5PRIME_UTR','3PRIME_UTR','NMD_TRANSCRIPT','SPLICE_SITE','STOP_GAINED','ESSENTIAL_SPLICE_SITE','WITHIN_MATURE_miRNA','STOP_LOST','PARTIAL_CODON','CODING_UNKNOWN']
#tags=['WITHIN_NON_CODING_GENE','INTRONIC','DOWNSTREAM','UPSTREAM','INTERGENIC','REGULATORY_REGION','NON_SYNONYMOUS_CODING','5PRIME_UTR','SYNONYMOUS_CODING','3PRIME_UTR','NMD_TRANSCRIPT','ESSENTIAL_SPLICE_SITE','SPLICE_SITE','STOP_GAINED','STOP_LOST','WITHIN_MATURE_miRNA']
tags=['regulatory_region_variant','downstream_gene_variant','upstream_gene_variant','intron_variant','nc_transcript_variant','non_coding_exon_variant','intergenic_variant','feature_elongation','feature_truncation','TF_binding_site_variant','splice_region_variant','missense_variant','synonymous_variant','3_prime_UTR_variant','NMD_transcript_variant','5_prime_UTR_variant','splice_donor_variant','inframe_deletion','splice_acceptor_variant','frameshift_variant','initiator_codon_variant','inframe_insertion','stop_lost','TFBS_ablation','stop_gained','stop_retained_variant','coding_sequence_variant','mature_miRNA_variant','incomplete_terminal_codon_variant','transcript_ablation']
diktag_snp={}
diktag_indel={}
for cat in tags: 
	diktag_snp[cat]=0
	diktag_indel[cat]=0

#this work by chr!!
print chr 
#pat tho vcf 
inputfile=gzip.open('/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1/chr_split/fvg.vqsr.beagle.impute2.anno.20140109.csq.pop.vcf.gz.%s.vcf.gz' %(chr), 'r')

for line in inputfile:
	x=line.split()
	if chr==23:
		chr="X"
 
	if x[0]==chr:
		#print line  
		for cat in tags:
			if re.search(cat, line): 
				if re.search('INDEL', line): diktag_indel[cat]+=1
				else: diktag_snp[cat]+=1
out=open('/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1/chr_split/FUNKANN/funkann.chr%s.table' %(chr), 'w')
sys.stdout=out
print 'CAT\tINDEL\tSNP' 
					
for cat in tags: print '%s\t%s\t%s' %(cat,diktag_indel[cat], diktag_snp[cat])

