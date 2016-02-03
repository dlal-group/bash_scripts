#!/usr/local/bin/bash
#using a job array approach 

# Fix multiallelic issue, trimming all sites with alternative allele uncounted
set -e
region=${1}
file=${2}
filename=`basename ${file}`
basedir=`dirname ${file}`

#now use each line of the chunk file to trim the file by region
bcftools view -a -r ${region} -O z -o ${basedir}/TRIMMED/M3.${filename}.${region}.vcf.gz ${file}
tabix -p vcf ${basedir}/TRIMMED/M3.${filename}.${region}.vcf.gz

