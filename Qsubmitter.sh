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
tar -xzvf MERGER.tgz MERGER/ALL/CHR${chr}/chr${chr}.geno.gz

