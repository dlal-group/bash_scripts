#!/usr/local/bin/bash
#reapply recalibration with a selected threshold/filter passed sites in a chr aware way
# Args
in_file=$1 #/lustre/scratch113/projects/fvg_seq/variant_calling/20130828/pooled/22.vcf.gz
out_path=$2 #/lustre/scratch113/projects/fvg_seq/variant_calling/20130828/pooled/22.vqsr.vcf.gz.part.vcf.gz
sens_threshold=$3 #this should be the minVQSlod value selected 
# 10.0453 => 98% sensitivity treshold
# 8.6247 => 99% sensitivity treshold

mode=$4

out_file=`basename ${in_file}`
# chr=`echo ${out_file%%.*}`

#let's use bcftools to annotate and change the annotation on vqsrfilter
# /nfs/users/nfs_m/mc14/Work/bash_scripts/ja_runner_par.sh /nfs/users/nfs_m/mc14/Work/bash_scripts/vqrs_reapply.sh vcf_files.list /nfs/users/nfs_m/mc14/fvg_seq/variant_refinemet/20140213_VQSR_REFILTER 8.6247 SNP
# bcftools2 filter /lustre/scratch113/projects/fvg_seq/variant_calling/20130828/pooled/10.vqsr.vcf.gz -m + -i INFO/VQSLOD>8.6247 -O v | bgzip -c  > /lustre/scratch113/projects/fvg_seq/variant_calling/20130828/pooled/10.vqsr.vcf.gz

# bcftools filter -i '%FILTER="PASS"' -O v ${in_file} | bgzip -c  > ${out_path}/${out_file}.filt.vcf.gz
bcftools filter -m + -i "INFO/VQSLOD>${sens_threshold}" ${in_file} -O v -o ${out_path}/${out_file}.filt.vcf.gz

#index the file
tabix -p vcf ${out_path}/${out_file}.filt.vcf.gz