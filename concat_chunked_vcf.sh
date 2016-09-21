#!/usr/local/bin/bash

#concat vcf files back together
set -e
file=${1}
filename=`basename ${file}`
basedir=`dirname ${file}`

# /nfs/users/nfs_m/mc14/carl_seq/variant_refinement/20092015_IMPUTE2/TRIMMED/M3.tt.1.vcf.gz.*.vcf.gz
ls ${basedir}/TRIMMED/M3.tt.${filename}.*.vcf.gz > ${basedir}/TRIMMED/M3.${filename}.chunked_files

#now use each line of the chunk file to trim the file by region
bcftools concat -a -f ${basedir}/TRIMMED/M3.${filename}.chunked_files -O z -o ${basedir}/TRIMMED/M3.${filename}
tabix -p vcf ${basedir}/TRIMMED/M3.${filename}

for line in $(cat ${basedir}/TRIMMED/M3.${filename}.chunked_files)
do
	rm ${line}
	rm ${line}.tbi
done

