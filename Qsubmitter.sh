#this function needs some arguments:
# $1=Job name
# $2=out_file
# $3=error_file
# $4=working directory
# $5=command to execute

#!/bin/bash
#$ -S /bin/bash
#$ -N "extract_chr${chr}"
#$ -o "extract_chr${chr}.o"
#$ -e "extract_chr${chr}.e"
#$ -cwd
#$ -q all.q

chr=$1
#10/02/2016
#extract file from a tar archive
# tar -xzvf MERGER.tgz MERGER/ALL/CHR${chr}/chr${chr}.geno.gz

#generate map file
cut -f 2,3,5 -d " " chr${chr}.geno_info > chr${chr}.part1_map

zcat chr${chr}.geno.gz | cut -f 2-5 -d " " > chr${chr}.part2_map

awk 'FNR==NR{a[$1,$2]=$3;next}{if(a[$1,$2]) print $0,a[$1,$2];else print $0,"NA"}' chr${chr}.part1_map chr${chr}.part2_map > chr${chr}.map