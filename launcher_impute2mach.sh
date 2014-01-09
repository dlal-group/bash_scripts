#!/usr/local/bin/bash

#launcher for a job array for impute2mach conversion

if [ $# -lt 4 ]
then
        echo "Missing arguments!!"
        echo "Usage:"
        echo "launcher_impute2mach.sh <chr> <geno_info_file_path> <sample_file> <output_file_path>"
        exit 1
fi

chr=$1
geno_file_path=$2
info_file_path=$2
sample_file=$3
output_file_path=$4

#Args submitted
echo ${chr}
echo ${geno_file_path}
echo ${info_file_path}
echo ${sample_file}
echo ${output_file_path}


geno_files=(`ls $geno_file_path/*.gen.gz`);
info_files=(`ls $info_file_path/*.gen_info`);
length=${#geno_files[@]};
bsub -J "convert[1-$length]" -o "%J_convert.log" -e "%J_convert.err" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q basement -- /nfs/users/nfs_m/mc14/Work/bash_scripts/impute2mach_launcher.sh $chr $geno_file_path $info_file_path $sample_file $output_file_path
