#!/usr/local/bin/bash
vcfchk_path=$1

(echo -e "CHR\tVARIANTS\tSNPS\tSingleton_SNPS\tMONO_SNPS\tINDELS\tSingleton_INDELS\tMONO_INDELS\tMULTI_SITES\tMULTI_SNPS\tMULTI_INDELS";
for i in {1..22}
do
snps=`grep "^SN" ${vcfchk_path}/${i}.vcf.gz.vchk| grep "number of SNPs"|cut -f 4`
indels=`grep "^SN" ${vcfchk_path}/${i}.vcf.gz.vchk| grep "number of indels"|cut -f 4`
sis_snps=`grep "^SiS" ${vcfchk_path}/${i}.vcf.gz.vchk | cut -f 4`
sis_indels=`grep "^SiS" ${vcfchk_path}/${i}.vcf.gz.vchk| cut -f 7`;
multi_sites=`grep "^SN" ${vcfchk_path}/${i}.vcf.gz.vchk| grep "number of multiallelic sites"|cut -f 4`;
multi_snps=`grep "^SN" ${vcfchk_path}/${i}.vcf.gz.vchk| grep "number of multiallelic SNP sites"|cut -f 4`;
multi_indel=$[multi_sites - multi_snps]
tot_v=$[snps + indels]
mono_snp=`grep "^AF" ${vcfchk_path}/${i}.vcf.gz.vchk| awk '$3==0' |cut -f 4`
mono_indels=`grep "^AF" ${vcfchk_path}/${i}.vcf.gz.vchk| awk '$3==0' |cut -f 7`
tot_mono=$[mono_snp + mono_indels]
echo -e "${i}\t${tot_v}\t${snps}\t${sis_snps}\t${mono_snp}\t${indels}\t${sis_indels}\t${mono_indels}\t${multi_sites}\t${multi_snps}\t${multi_indel}";
done)
