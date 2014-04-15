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
chr=`echo ${out_file%%.*}`

# /software/jre1.7.0_25/bin/java -Xmx12000m -Xms12000m -Xss280m -server -XX:+UseSerialGC -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-2.7-2-g6bda569/GenomeAnalysisTK.jar \
# -T ApplyRecalibration -U LENIENT_VCF_PROCESSING --ts_filter_level ${sens_threshold} --mode ${mode} \
# -R /lustre/scratch111/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa \
# --input ${in_file} \
# -recalFile /lustre/scratch113/projects/fvg_seq/variant_calling/20130828/pooled/vqsr.sites.snps.vcf.gz \
# -tranchesFile /lustre/scratch113/projects/fvg_seq/variant_calling/20130828/pooled/vqsr.sites.snps.tranches \
# -o ${out_path}/${out_file}.vqsr.vcf
#let's use bcftools to annotate an change the annotation on vqsrfilter
# /nfs/users/nfs_m/mc14/Work/bash_scripts/ja_runner_par.sh /nfs/users/nfs_m/mc14/Work/bash_scripts/vqrs_reapply.sh vcf_files.list /nfs/users/nfs_m/mc14/fvg_seq/variant_refinemet/20140213_VQSR_REFILTER 8.6247 SNP
# bcftools2 filter /lustre/scratch113/projects/fvg_seq/variant_calling/20130828/pooled/10.vqsr.vcf.gz -m + -i INFO/VQSLOD>8.6247 -O v | bgzip -c  > /lustre/scratch113/projects/fvg_seq/variant_calling/20130828/pooled/10.vqsr.vcf.gz

# bcftools2 filter -m + -i "INFO/VQSLOD>${sens_threshold}" -O v ${in_file} | bgzip -c  > ${out_path}/${out_file}
bcftools2 filter -i '%FILTER="PASS"' -O v ${in_file} | bgzip -c  > ${out_path}/${out_file}.filt.vcf.gz

#index the file
tabix -p vcf ${out_path}/${out_file}.filt.vcf.gz