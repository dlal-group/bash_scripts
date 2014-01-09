#!/usr/local/bin/bash

#commands to extract list of sites and other stuff needed for the gemma_additional_calc scrips
#args:
#$1 : trait
#$2 : result path
#$3 : outpath

#works by chr
trait=$1
result_path=$2
outpath=$3

mkdir -p ${outpath}

if [ $LSB_JOBINDEX -eq 23 ]
then
        chr="X"
else
        chr=$LSB_JOBINDEX
fi

#extract a list of snps to pass to snptest for stat calc
fgrep -v "beta" ${result_path}/${trait}.chr${chr}.tab.assoc.txt | cut -f 2|tr "-" " " |awk '{if($2=="NA") print $1;else print $2}' > ${outpath}/${trait}.chr${chr}.rs_to_include

