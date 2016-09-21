#!/bin/bash
#$ -S /bin/bash
#$ -N "extract_chr"
#$ -o "$JOB_ID_extract_chr.o"
#$ -e "$JOB_ID_extract_chr.e"
#$ -cwd
#$ -q all.q

set -e
chr=$1
#10/02/2016
#extract file from a tar archive
# tar -xzvf MERGER.tgz MERGER/ALL/CHR${chr}/chr${chr}.geno.gz

#generate map file
# cut -f 2,3,5 -d " " /netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.geno_info > /netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.part1_map
# zcat /netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.geno.gz | cut -f 2-5 -d " " > /netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.part2_map
# (echo "SNP Position A0 A1 Rsq";awk 'FNR==NR{a[$1,$2]=$3;next}{if(a[$1,$2]) print $0,a[$1,$2]}' /netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.part1_map /netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.part2_map) | tr " " "\t"> /netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.map

#17/02/2016
#extract info for results qc after GWA
# out_path=$2
# cohort=$3
# out_suffix=$4
# # mkdir -p ${out_path}/${chr}/
# # tar -xzvf ${cohort}_MetS_score_${chr}_*.tar.gz -C ${out_path}/${chr}/


# R CMD BATCH '--args '${cohort}' '${chr}' '${out_path}' '${out_path}/${chr}/${cohort}_MetS_score_${chr}_${out_suffix}.csv'' ~/scripts/r_scripts/GWA_qc.r ${out_path}/${chr}/${cohort}_MetS_score_analysis_chr${chr}.Rout

#17/02/2016
# extract certainty from info files and add it to GWA results
file1=/netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.geno_info
file2=/home/cocca/analyses/MetabolicSyndrome/VBI/results/${chr}/VBI_MetS_score_${chr}_Feb_16_2016_cocca_results.csv

awk 'BEGIN{OFS=","}{print $2,$3,$6}' ${file1} > ${file1}.cert.txt
#join results and certainty info
header1=`head -1 ${file2}`
(echo "${header1},avpostprob";awk 'BEGIN{FS=",";OFS=","}FNR==NR{a[$1,$2]=$3;next}{if(a[$1,$3]) print $0,a[$1,$3]}' ${file1}.cert.txt ${file2}) > ${file2}.complete.csv
