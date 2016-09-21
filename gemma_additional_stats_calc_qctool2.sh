#!/usr/local/bin/bash

#Pipelene to extract some useful stats after a gemma analysis
#We want a file with :
#CHROM POS RSID A0 A1 MINOR_ALL MAJOR_ALL MAF QUALITY INFO
#for each dataset used in the analysis.
#we use qctool and info files from imputation to obtain the file

#define args
#$1=imputation file path 
#$2=sample file path
#$3=rsID list file path for snps exclusion
#$4=sample list file path for exclusion

imp_file_path=$1
sample_file_path=$2
rs_to_remove_path=$3
sample_to_remove_path=$4


sample_file_name=`basename ${sample_to_remove_path}`
echo $sample_file_name
stat_filename=${sample_file_name%%_*}
echo $stat_filename
trait=${stat_filename#*.}
echo $trait
out_path=$5/$trait
echo $out_path
mkdir -p $out_path/LOGS

for chr in {1..22}
#for chr in $6
do

#First step:
#Samples exclusion

#bsub -J "s1_${trait}_chr${chr}" -o "$out_path/LOGS/%J_s1_chr${chr}.log" -e "$out_path/LOGS/%J_s1_chr${chr}.err" -G "team151" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q basement qctool -g $imp_file_path/CHR${chr}/chr${chr}.gen.gz -s $sample_file_path -excl-samples $sample_to_remove_path -og $out_path/chr${chr}_sample_filtered.gen.gz -os $out_path/Filtered_chr${chr}_samples.sample
bsub -J "s1_${trait}_chr${chr}" -o "$out_path/LOGS/%J_s1_chr${chr}.log" -e "$out_path/LOGS/%J_s1_chr${chr}.err" -G "team151" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q basement qctool -g $imp_file_path/CHR${chr}/chr${chr}.gen.gz -s $sample_file_path -assume-chromosome ${chr} -excl-samples $sample_to_remove_path -og $out_path/chr${chr}_sample_filtered.bgen -os $out_path/Filtered_chr${chr}_samples.sample

#Second step:
#SNPs exclusion:exclude all snps not used in the analyses and stats calc

#bsub -J "s2_${trait}_chr${chr}" -w "ended(s1_${trait}_chr${chr})" -o "$out_path/LOGS/%J_s2_chr${chr}.log" -e "$out_path/LOGS/%J_s2_chr${chr}.err" -G "team151" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q basement qctool -g $out_path/chr${chr}_sample_filtered.gen.gz -s $out_path/Filtered_chr${chr}_samples.sample -incl-rsids $rs_to_remove_path -os $out_path/filtered_${chr}_samples.stats.sample -snp-stats $out_path/${stat_filename}.${chr}.stats
bsub -J "s2_${trait}_chr${chr}" -w "ended(s1_${trait}_chr${chr})" -o "$out_path/LOGS/%J_s2_chr${chr}.log" -e "$out_path/LOGS/%J_s2_chr${chr}.err" -G "team151" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q basement qctool -g $out_path/chr${chr}_sample_filtered.bgen -s $out_path/Filtered_chr${chr}_samples.sample -assume-chromosome ${chr} -incl-rsids $rs_to_remove_path -os $out_path/filtered_${chr}_samples.stats.sample -snp-stats $out_path/${stat_filename}.${chr}.stats

#Third step:
#Stats extraction

#bsub -J "s3_${trait}_chr${chr}" -w "ended(s2_${trait}_chr${chr})" -o "$out_path/LOGS/%J_s3_chr${chr}.log" -e "$out_path/LOGS/%J_s3_chr${chr}.err" -G "team151" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q basement qctool -g $out_path/chr${chr}_sample_filtered.gen.gz -s $out_path/Filtered_chr${chr}_samples.sample -assume-chromosome ${chr} -snp-stats $out_path/${stat_filename}.${chr}.stats -os $out_path/filtered_${chr}_samples.stats.sample

#now create the output file we want:
#first for each chromosome, then we merge all toghether
#We want a file with :
#CHROM POS RSID A0 A1 MINOR_ALL MAJOR_ALL MAF QUALITY INFO-

#reorder and pick only the right columns
bsub -J "s4_${trait}_chr${chr}" -w "ended(s2_${trait}_chr${chr})" -o "$out_path/LOGS/%J_s4_chr${chr}.log" -e "$out_path/LOGS/%J_s4_chr${chr}.err" -G "team151" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q basement /nfs/users/nfs_m/mc14/Work/bash_scripts/reorder_gemma_stats.sh $out_path ${stat_filename}.${chr}.stats $imp_file_path/CHR${chr}/chr${chr}.gen_info

#once finished step 2 we can make a little cleanup and compression
#bsub -J "s5_${trait}_chr${chr}" -w "ended(s2_${trait}_chr${chr})" -o "$out_path/LOGS/%J_s5_chr${chr}.log" -e "$out_path/LOGS/%J_s5_chr${chr}.err" -G "team151" -q normal "rm $out_path/chr${chr}_snp_filtered.gen.gz"

#bsub -J "s6_${trait}_chr${chr}" -w "ended(s2_${trait}_chr${chr})" -o "$out_path/LOGS/%J_s6_chr${chr}.log" -e "$out_path/LOGS/%J_s6_chr${chr}.err" -G "team151" -q normal "rm $out_path/chr${chr}_sample_filtered.gen.gz"
bsub -J "s6_${trait}_chr${chr}" -w "ended(s2_${trait}_chr${chr})" -o "$out_path/LOGS/%J_s6_chr${chr}.log" -e "$out_path/LOGS/%J_s6_chr${chr}.err" -G "team151" -q normal "rm $out_path/chr${chr}_sample_filtered.bgen"

#once finished step 4 we can make a little cleanup and compression
bsub -J "s5_${trait}_chr${chr}" -w "ended(s4_${trait}_chr${chr})" -o "$out_path/LOGS/%J_s5_chr${chr}.log" -e "$out_path/LOGS/%J_s5_chr${chr}.err" -G "team151" -q normal "rm $out_path/${stat_filename}.${chr}.stats;rm $out_path/Filtered_chr${chr}_samples.sample"

done
