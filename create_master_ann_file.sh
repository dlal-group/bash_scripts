#!/usr/local/bin/bash

# This is the runner file run by bsub
# Arguments: 
#filepath
#chr
#outpath

# Environment variables: LSB_JOBINDEX
# index=$[LSB_JOBINDEX - 1]

# Modified to be run as a job array
# files=(`ls $1/chr${2}.*.gen.gz`)

# file=${files[$index]}
file=$1

sample_f=${file%.*}_samples

# added to provide right format for sample files
sed -i 's/-9/NA/g' ${sample_f}

# out_path=$3
out_path=$2

mkdir -p ${out_path}

out_name=`basename ${file}`

# /nfs/team151/software/bin/snptest_v2 -summary_stats_only -data ${file} ${sample_f} -o ${out_path}/${out_name}.stats -hwe -log ${out_path}/${out_name}.log 
#Use the new version to work on chrX
# /nfs/team151/software/snptest_v2.5-beta4_Linux_x86_64_static/snptest_v2.5-beta4 -summary_stats_only -data ${file} ${sample_f} -assume_chromosome $2 -o ${out_path}/${out_name}.stats -hwe -log ${out_path}/${out_name}.log 
snptest_v2 -summary_stats_only -data ${file} ${sample_f} -o ${out_path}/${out_name}.stats -assume_chromosome X -hwe -log ${out_path}/${out_name}.log 

gzip ${out_path}/${out_name}.stats