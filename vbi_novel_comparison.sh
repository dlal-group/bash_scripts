#!/usr/local/bin/bash

#script to check what in novel vbi sites is also in uk10k set
#what we need:
#$1=input_1_path
#$2=input_2_path
#$3=chr number
#$4=output_path

if [ $5 != 'overlap' ]
then
	#EXTRACT NOT OVERLAPPING LIST
	mkdir -p $4/NOT_OVERLAP
	#cat <(zcat $2/UK10K.chr$3.snps.vcf.gz | grep ^# ) <(fgrep -w -v -f <(cut -f 1,2 $1/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$3.re_ann.NOT_OVERLAP.NO_RSID.map ) <(tabix $2/UK10K.chr$3.snps.vcf.gz chr $3)) | bgzip -c > $4/NOT_OVERLAP/UK10K.chr$3.snps.not_overlap.vcf.gz 

	#Modified 27/10/2012 in order to remake all checks
	cat <(zcat $2/UK10K.chr$3.snps.vcf.gz | grep ^# ) <(fgrep -w -v -f <(cut -f 1,2 $1/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$3.re_ann.NOT_OVERLAP.NO_RSID.former_SNPs_removed.map ) <(tabix $2/UK10K.chr$3.snps.vcf.gz chr $3)) | bgzip -c > $4/NOT_OVERLAP/UK10K.chr$3.snps.not_overlap.vcf.gz 
	tabix $4/NOT_OVERLAP/UK10K.chr$3.snps.not_overlap.vcf.gz
else
	#EXTRACT OVERLAPPING LIST
	mkdir -p $4/OVERLAP
	#cat <(zcat $2/UK10K.chr$3.snps.vcf.gz | grep ^# ) <(fgrep -w -f <(cut -f 1,2 $1/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$3.re_ann.NOT_OVERLAP.NO_RSID.map ) <(tabix $2/UK10K.chr$3.snps.vcf.gz chr $3)) | bgzip -c > $4/OVERLAP/UK10K.chr$3.snps.overlap.vcf.gz

	#Modified 27/10/2012 in order to remake all checks
	cat <(zcat $2/UK10K.chr$3.snps.vcf.gz | grep ^# ) <(fgrep -w -f <(cut -f 1,2 $1/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$3.re_ann.NOT_OVERLAP.NO_RSID.former_SNPs_removed.map ) <(tabix $2/UK10K.chr$3.snps.vcf.gz chr $3)) | bgzip -c > $4/OVERLAP/UK10K.chr$3.snps.overlap.vcf.gz
	tabix $4/OVERLAP/UK10K.chr$3.snps.overlap.vcf.gz
fi
