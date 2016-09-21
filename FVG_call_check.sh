#!/usr/local/bin/bash

#Pipeline script to performs check on vcf files
#
#Steps involved:
#Extract genotypes from plink files
chr=$1
inpath=$2
outdir=$3
keeplist=$4
mkdir -p ${outdir}

# plink --noweb --bfile /lustre/scratch113/projects/uk10k/users/jh21/imputed/fvg/fvg_370/shapeit/chr5 --chr 5 --keep ../test.keeplist --make-bed --out test_chr5_single_sample
plink --noweb --bfile ${inpath}/chr${chr} --chr ${chr} --keep ${keeplist} --make-bed --out ${outdir}/FVG_seq_subset_chr${chr}
