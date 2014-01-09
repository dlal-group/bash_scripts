#!/usr/local/bin/bash

chr=${LSB_JOBINDEX}
path=$1

#fgrep -w -f <( cut -f 2 ~/Work/SANGER/DOCS/toploci_WGS.txt ) <( tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.vcf.gz ${chr}) > looked_chr${chr}.txt
fgrep -w -f <( cut -f 2 ~/Work/SANGER/DOCS/toploci_WGS.txt ) $path/CHR${chr}/chr${chr}.gen_info_rs_only_maf > looked_chr${chr}.txt
