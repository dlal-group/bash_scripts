#!/usr/local/bin/bash

#Pipeline script to performs check on vcf files
#
#Steps involved:
#Extract genotypes from plink files
inpath=$1
outdir=$2
keeplist=$3
mkdir -p ${outdir}

for chr in {1..22} X
do
	bsub -J "plink_run_${chr}" -o "${outdir}/%J_plink_run_${chr}.o" -M 3000 -R"select[mem>3000] rusage[mem=3000]" -q normal -- plink --noweb \
	--bfile ${inpath}/chr${chr} \
	--chr ${chr} \
	--keep ${keeplist} \
	--make-bed \
	--out ${outdir}/FVG_seq_subset_chr${chr}
done
