#!/usr/local/bin/bash

#script to launch impute2bimbam.sh 
#args:
#$1=chromosome files geno path
#$2=transformed file output

if [ $# -lt 2 ]
then
        echo "ERROR!!Missing arguments!!"
        echo "USAGE:"
        echo "impute2bimbam_launcher.sh < full chr file path> <full output file path>"
exit 1
fi


# for chr in {1..22} X
for chr in 15
do

# bsub -J "impute2bimbam_${chr}" -o "%J_impute2bimbam_${chr}.log" -M1000 -R"select[mem>1000] rusage[mem=1000]" \
# -q normal impute2bimbam.sh ${chr} $1 $2
impute2bimbam.sh ${chr} $1 $2

done