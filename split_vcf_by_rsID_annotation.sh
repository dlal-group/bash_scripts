#!/usr/local/bin/bash
#
#Script to split vcf files in missing rsID annotation vcf

seq_file_path=$2
splitted_seq_file_path=$3

mkdir -p RS_ann_splitted/NO_RSID

#for i in {1..22}
#do
zcat ${splitted_seq_file_path}/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.vcf.gz | head -200 | grep '^#' > header_chr$1.tmp
tabix -f ${splitted_seq_file_path}/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.vcf.gz
tabix ${splitted_seq_file_path}/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.vcf.gz chr $1 | awk '{OFS="\t"}{if ($3 == ".") print $0}' > splitted_ann_chr$1.tmp
cat header_chr$1.tmp splitted_ann_chr$1.tmp | bgzip -c > RS_ann_splitted/NO_RSID/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.no_rsIDs.vcf.gz

rm header_chr$1.tmp
rm splitted_ann_chr$1.tmp
# done
