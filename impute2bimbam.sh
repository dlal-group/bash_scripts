#!/usr/local/bin/bash

#script to convert impute2 output files to bimbam comma sepatarated mean genotype files
#args:
#$1=chr number
#$2=chromosome files geno path
#$3=transformed file output

#convert first each chunk
# files=`ls $2/chr$1.*.gen.gz`

# for file in $files
# do
# mkdir -p $3
# filename=`basename $file`
# bsub -J "impute2bimbam_$filename" -o "%J_impute2bimbam_$filename.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
# -q basement impute2bimbam_inner.sh $file $filename $3

# done

#now launch conversion for the huge whole chr file
chr=$1
whole_file=`ls $2/chr$1.gen.gz`
mkdir -p $3
wholefile_name=`basename $whole_file`
# bsub -J "impute2bimbam_$wholefile_name" -o "%J_impute2bimbam_$wholefile_name.log" -M16000000 -R"select[mem>16000] rusage[mem=16000]" -q basement impute2bimbam_inner.sh $whole_file $wholefile_name $3
bsub -J "impute2bimbam_$wholefile_name" -o "$3/%J_impute2bimbam_$wholefile_name.o" -M7000 -R"select[mem>7000] rusage[mem=7000]" -q normal impute2bimbam_inner.sh $whole_file ${chr} $3

