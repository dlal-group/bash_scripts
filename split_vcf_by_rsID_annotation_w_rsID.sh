#!/usr/local/bin/bash
#
#Script to split vcf files in rsId annotated vcf and missing rsID annotation vcf

seq_file_path=$2
splitted_seq_file_path=$3

#mkdir -p RS_ann_splitted

#for i in {1..22}
#do
zcat ${seq_file_path}/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.vcf.gz | head -200 | grep '^#' > header_w_rs_chr$1.tmp
tabix -f ${splitted_seq_file_path}/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.vcf.gz
tabix ${splitted_seq_file_path}/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.vcf.gz chr $1 | awk '{OFS="\t"}{if ($3 != ".") print $0}' > splitted_w_rs_ann_chr$1.tmp
cat header_w_rs_chr$1.tmp splitted_w_rs_ann_chr$1.tmp | bgzip -c > RS_ann_splitted/W_RSID/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.w_rsIDs.vcf.gz

rm header_w_rs_chr$1.tmp
rm splitted_w_rs_ann_chr$1.tmp
# done
