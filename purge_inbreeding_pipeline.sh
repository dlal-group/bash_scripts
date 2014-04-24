#!/usr/local/bin/bash
#Script to create a pipeline for the purge through inbreeding work
#
#

# Merge different popuplation together

#Extract IBD information for each population. Here we are using plink2 (1.9)
pops="FVG VBI TSI CEU"
pops="VBI TSI CEU"
for pop in $pops
do

case $pop in
  FVG )
    pop_path=/lustre/scratch113/projects/fvg_seq/20140410/INGI/FVG
    ;;
  VBI )
    pop_path=/lustre/scratch113/projects/fvg_seq/20140410/INGI/VBI
      ;;
  TSI )
    pop_path=/lustre/scratch113/projects/fvg_seq/20140410/TGP/TSI
      ;;
  CEU )
    pop_path=/lustre/scratch113/projects/fvg_seq/20140410/TGP/CEU
      ;;
esac
  bsub -J"ibd_${pop}" -o"%J_ibd_${pop}.o" -q yesterday -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --vcf ${pop_path}/22.vcf.gz --double-id --biallelic-only --genome gz --parallel 1 2 --out ibd_${pop}
done
