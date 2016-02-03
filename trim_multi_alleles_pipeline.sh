#!/usr/local/bin/bash
#using a job array approach 

# Fix multiallelic issue, trimming all sites with alternative allele uncounted
set -e
file=$1
filename=`basename ${file}`
basedir=`dirname ${file}`

mkdir -p ${basedir}/TRIMMED/

#we need to select the multialleic sites, first, since we just want to work on them
bcftools view -m 3 ${file} -O z -o ${basedir}/TRIMMED/M3.tt.${filename}
tabix -p vcf -f ${basedir}/TRIMMED/M3.tt.${filename}

#now create chunks for each chr
python ~/Work/bash_scripts/create_chunks.py ${basedir}/TRIMMED/M3.tt.${filename}

#now use each line of the chunk file to trim the file by region
mkdir -p ${basedir}/TRIMMED/LOGS;size=`wc -l ${basedir}/TRIMMED/M3.tt.${filename}.chunks |cut -f 1 -d " "`; bsub -J "trim_${filename}[1-${size}]" -o "${basedir}/TRIMMED/LOGS/%J_trim_${filename}.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/trim_multi_alleles.sh ${basedir}/TRIMMED/M3.tt.${filename}.chunks ${basedir}/TRIMMED/M3.tt.${filename}
#${basedir}/TRIMMED/M3.tt.1.vcf.gz

#in the meanwhile we create a file for all the remaining sites: we're not goin to work on them!
bsub -J "extract_M2_${filename}" -o "${basedir}/TRIMMED/LOGS/%J_extract_M2_${file_name}.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- bcftools view -M 2 -O z -o ${basedir}/TRIMMED/M2.${filename} ${file}
bsub -J "index_M2_${filename}" -o "${basedir}/TRIMMED/LOGS/%J_index_M2_${file_name}.o" w "ended(extract_M2_${filename})" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- tabix -p vcf -f ${basedir}/TRIMMED/M2.${filename}

#concat with the -a options, sort stuff out by itself!!
#We need to wait for all the jobs to end, than we concat M3 chunks back together
bsub -J "concat_M3_${filename}" -o "${basedir}/TRIMMED/LOGS/%J_concat_M3_${file_name}.o" w "ended(trim_${filename}*)" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/concat_chunked_vcf.sh ${file}

#We need to wait for all the jobs to end, than we concat M2 and M3 stuff back together
bsub -J "concat_M2M3_${filename}" -o "${basedir}/TRIMMED/LOGS/%J_concat_M2M3_${file_name}.o" w "ended(concat_M3_${filename})" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- bcftools concat -a ${basedir}/TRIMMED/M2.${filename} ${basedir}/TRIMMED/M3.${filename} -O z -o ${basedir}/TRIMMED/${filename}
bsub -J "index_M2M3_${filename}" -o "${basedir}/TRIMMED/LOGS/%J_index_M2M3_${file_name}.o" w "ended(concat_M2M3_${filename})" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- tabix -p vcf -f ${basedir}/TRIMMED/${filename}

echo "rm ${basedir}/TRIMMED/M3.tt.${filename};rm ${basedir}/TRIMMED/M3.tt.${filename}.tbi;rm ${basedir}/TRIMMED/M2.${filename};rm ${basedir}/TRIMMED/M3.${filename};rm ${basedir}/TRIMMED/M2.${filename}.tbi;rm ${basedir}/TRIMMED/M3.${filename}.tbi" | bsub -J "clean_${filename}" -o "${basedir}/TRIMMED/LOGS/%J_clean_${file_name}.o" w "ended(index_M2M3_${filename})" -M1000 -R"select[mem>1000] rusage[mem=1000]" -q normal
