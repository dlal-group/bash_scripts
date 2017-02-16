#!/usr/bin/env bash
set -e

if [ $# -lt 5 ]
then
	echo "Missin' arguments!!"
	echo "Usage:"
	echo "impute2mach_launcher.sh <chr> <geno_file> <info_file> <sample_file> <output_file>"
	exit 1
fi

chr=$1
geno_file=$2
info_file=$3
sample_file=$4
output_file_path=$5

# index=$[LSB_JOBINDEX - 1]

#Args submitted
echo ${chr}
echo ${geno_file}
echo ${info_file}
echo ${sample_file}
echo ${output_file_path}

# geno_files=(`ls $geno_file_path/*.gen.gz`);
# info_files=(`ls $info_file_path/*.gen_info`);

#we need to be in the same folder of the output if we want to avoid errors
mkdir -p ${output_file_path}/${chr}
cd ${output_file_path}/${chr}/

# geno_file=${geno_files[$index]}
# info_file=${info_files[$index]}
output_file=${output_file_path}/${chr}

#Extract geno file
geno_file_name=`basename $geno_file`
#11/10/2016 - skip duplicate romoval: we didit already in the previous step
# gunzip -c ${geno_file} | awk '!($0 in a){a[$0];print}' > ${output_file_path}/${chr}/${geno_file_name}.gen
gunzip -c ${geno_file}  > ${output_file_path}/${chr}/${geno_file_name}.gen


echo "Convert to file vector..."
#script to launch with qsub impute2mach conversion
R CMD BATCH '--args '${chr}' '${geno_file_name}.gen' '${info_file}' '${sample_file}' '${output_file}'' /home/cocca/scripts/r_scripts/impute2mach.R ${output_file_path}/${chr}.Rout

echo "Conversion ended!"
# echo "Sort the converted file..."
#sort the files obtained in mach format
output_file_suffix_path=${output_file}

#clean useless files
rm -r ${output_file_path}/${chr}/
# sort -k1 ${output_file_suffix_path}.machdose -o ${output_file_suffix_path}.machdose 

#add the 1-> field in the dose file
# awk '{OFS=" "}{print "1->"$0}' ${output_file_suffix_path}.machdose > ${output_file_suffix_path}_sorted.dose
