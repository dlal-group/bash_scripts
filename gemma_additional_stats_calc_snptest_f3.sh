#!/usr/local/bin/bash

#Pipelene to extract some useful stats after a gemma analysis
#We want a file with :
#CHROM POS RSID A0 A1 MINOR_ALL MAJOR_ALL MAF QUALITY INFO
#for each dataset used in the analysis.
#we use snptest and info files from imputation to obtain the file

#define args
#$1=imputation file path 
#$2=sample file path
#$3=rsID list file path for snps inclusion
#$4=sample list file path for exclusion
#$5=outh path

if [ $# -lt 5 ]
then
	echo "missing arguments!!"
	echo "USAGE: gemma_additional_stats_calc_snptest.sh <imputed_file_path> <sample_file_path> <rs_id_file_path> <sample_file_path> <out_path>"
	exit 1
fi

imp_file_path=$1
sample_file_path=$2
rs_to_include_path=$3
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

for chr in {1..22} X
#for chr in $6
do

#First step:
#Extract exclusion snps list
#Use as snps list the snps in the gemma output
#Snp list from gemma 

# bsub -J "s1_${trait}_chr${chr}" -o "$out_path/LOGS/%J_s1_chr${chr}.log" -e "$out_path/LOGS/%J_s1_chr${chr}.err" -G "team151" -M8000000 -R "select[mem>8000] rusage[mem=8000]" -q normal -- "fgrep -w -v -f $rs_to_include_path $imp_file_path/CHR${chr}/chr${chr}.gen_info | cut -f 2 -d ' ' | grep -v position > $out_path/rs_to_remove_${chr}.list"
bsub -J "s1_${trait}_chr${chr}" -o "$out_path/LOGS/%J_s1_chr${chr}.log" -e "$out_path/LOGS/%J_s1_chr${chr}.err" -G "team151" -M5000 -R "select[mem>5000] rusage[mem=5000]" -q normal -- "fgrep -w -v -f $rs_to_include_path $imp_file_path/chr${chr}.gen_info | cut -f 2 -d ' ' | grep -v position > $out_path/rs_to_remove_${chr}.list"

#Second step:
#Stats extraction

bsub -J "s2_${trait}_chr${chr}" -o "$out_path/LOGS/%J_s2_chr${chr}.log" -e "$out_path/LOGS/%J_s2_chr${chr}.err" -w "ended(s1_${trait}_chr${chr})" -G "team151" -M8000 -R "select[mem>8000] rusage[mem=8000]" -q normal -- /nfs/team151/software/bin/snptest_v2 -summary_stats_only -data $imp_file_path/chr${chr}.gen.gz $sample_file_path -exclude_snps $out_path/rs_to_remove_${chr}.list -exclude_samples $sample_to_remove_path -hwe -log $out_path/chr${chr}_filtered.log -o $out_path/${stat_filename}.${chr}.stats


#now create the output file we want:
#first for each chromosome, then we merge all toghether
#We want a file with :
#CHROM POS RSID A0 A1 MINOR_ALL MAJOR_ALL MAF QUALITY INFO-

#reorder and pick only the right columns
bsub -J "s3_${trait}_chr${chr}" -w "ended(s2_${trait}_chr${chr})" -o "$out_path/LOGS/%J_s3_chr${chr}.log" -e "$out_path/LOGS/%J_s3_chr${chr}.err" -G "team151" -M8000 -R "select[mem>8000] rusage[mem=8000]" -q normal -- /nfs/users/nfs_m/mc14/Work/bash_scripts/reorder_gemma_stats.sh $out_path $out_path/${stat_filename}.${chr}.stats $imp_file_path/chr${chr}.gen_info ${chr}

#once finished step 3 we can make a little cleanup and compression
#bsub -J "s4_${trait}_chr${chr}" -w "ended(s3_${trait}_chr${chr})" -o "$out_path/LOGS/%J_s4_chr${chr}.log" -e "$out_path/LOGS/%J_s4_chr${chr}.err" -G "team151" -q normal "rm $out_path/${stat_filename}.${chr}.stats"

#bsub -J "s6_${trait}_chr${chr}" -w "ended(s2_${trait}_chr${chr})" -o "$out_path/LOGS/%J_s6_chr${chr}.log" -e "$out_path/LOGS/%J_s6_chr${chr}.err" -G "team151" -q normal "rm $out_path/chr${chr}_sample_filtered.gen.gz"
#bsub -J "s6_${trait}_chr${chr}" -w "ended(s2_${trait}_chr${chr})" -o "$out_path/LOGS/%J_s6_chr${chr}.log" -e "$out_path/LOGS/%J_s6_chr${chr}.err" -G "team151" -q normal "rm $out_path/chr${chr}_sample_filtered.bgen"

#once finished step 4 we can make a little cleanup and compression
#bsub -J "s5_${trait}_chr${chr}" -w "ended(s4_${trait}_chr${chr})" -o "$out_path/LOGS/%J_s5_chr${chr}.log" -e "$out_path/LOGS/%J_s5_chr${chr}.err" -G "team151" -q normal "rm $out_path/${stat_filename}.${chr}.stats;rm $out_path/Filtered_chr${chr}_samples.sample"

done
