#!/bin/bash
#$ -S /bin/bash
#$ -N "extract_chr${chr}"
#$ -o "$JOB_ID_extract_chr${chr}.o"
#$ -e "$JOB_ID_extract_chr${chr}.e"
#$ -cwd
#$ -q all.q

set -e
chr=$1
#10/02/2016
#extract file from a tar archive
# tar -xzvf MERGER.tgz MERGER/ALL/CHR${chr}/chr${chr}.geno.gz

#generate map file
cut -f 2,3,5 -d " " /netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.geno_info > /netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.part1_map

zcat /netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.geno.gz | cut -f 2-5 -d " " > /netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.part2_map

(echo "SNP Position A0 A1 Rsq";awk 'FNR==NR{a[$1,$2]=$3;next}{if(a[$1,$2]) print $0,a[$1,$2]}' /netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.part1_map /netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.part2_map) | tr " " "\t"> /netapp/dati/daCinzia/Genotipi/VB_IMP_1KG/MERGER/ALL/CHR${chr}/chr${chr}.map