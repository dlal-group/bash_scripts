#!/usr/local/bin/bash
#reannotate vcf for multiallelic snps removal
#bsub -J "reannotate" -o "%J_reannotate.log" -M16000000 -R"select[mem>16000] rusage[mem=16000]" \
#-q basement "zcat esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz | vcf-subset -a | vcf-annotate --fill-ICF | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.reannotated.vcf.gz"

##update 08-15-2012
#this script was used to solve the rsIDs issue, reannotating the vcf files from vbi sequences
# set the base path for the annotation file (we assume is the current working dir)
annotation_base_path=`pwd`

if [ $# -eq 4 ]
then
	zcat $3/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.no_rsIDs.vcf.gz| vcf-annotate -a $annotation_base_path/dbSNP_splitted/chr$1_dbSNP_b137.tab.gz -c CHROM,POS,ID,-,-,-,- | bgzip -c > $2/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann_rsIDs.vcf.gz
else
	zcat $3/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.no_rsIDs.vcf.gz| vcf-annotate -a $annotation_base_path/dbSNP_splitted/chr$1_dbSNP_b137.tab.gz -c CHROM,POS,ID,REF,ALT,-,- | bgzip -c > $2/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann_rsIDs.vcf.gz

fi
