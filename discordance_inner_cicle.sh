#!/bin/bash

#conversion for seq files
plink --noweb --file $1 --chr $3 --recode --transpose --tab --allow-no-sex --out ${1}_chr$3
#conversion for other files
plink --noweb --file $2 --chr $3 --recode --transpose --tab --allow-no-sex --out ${2}_chr$3
#now launch the python script on any chr
#python2.7 ~/Work/bash_scripts/non_ref_discordance.py 
#370K/TEST_AD/biallelic_overl.SNP.unfilt.geno.seq.VB.chr20.tped
# SEQ/XX_sex_plink.overl_pos_NO_INDEL.esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.recode.vcf.chr20.tped
# ref_discordance_table.txt 
#20
python2.7 ~/Work/bash_scripts/non_ref_discordance.py \
${1}_chr$3.tped \
${2}_chr$3.tped \
$5 $3 $4
#move files
mkdir -p QC_OUT/CHR$3
mv *_chr$3.txt QC_OUT/CHR$3

