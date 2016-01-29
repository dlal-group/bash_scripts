#!/usr/local/bin/bash
#using a job array approach 

# Fix multiallelic issue, trimming all sites with alternative allele uncounted
set -e
filename=`basename ${file}`
basedir=`dirname ${file}`

mkdir -p ${basedir}/TRIMMED/

#we need to select the multialleic sites, first, since we just want to work on them
bcftools view -a -m 3 ${file} -O z -o ${basedir}/TRIMMED/M3.${filename}
tabix -p vcf ${basedir}/TRIMMED/M3.${filename}

#then we create a file for all the remaining sites: we're not goin to work on them!
bcftools view -M 2 -O z -o ${basedir}/TRIMMED/M2.${filename} ${file}
tabix -p vcf ${basedir}/TRIMMED/M2.${filename}

#concat with the -a options, sort stuff out by itself!!
bcftools concat -a ${basedir}/TRIMMED/M2.${filename} ${basedir}/TRIMMED/M3.${filename} -O z -o ${basedir}/TRIMMED/${filename}
tabix -p vcf ${basedir}/TRIMMED/${filename}

rm ${basedir}/TRIMMED/M2.${filename}
rm ${basedir}/TRIMMED/M3.${filename}
rm ${basedir}/TRIMMED/M2.${filename}.tbi
rm ${basedir}/TRIMMED/M3.${filename}.tbi
