#!/usr/local/bin/bash
#script to extract dosages for snps to use in the conditional analyses
#ARGS:
#$1 = genotypes file path
#$2 = snp list
#$3 = out folder 

#if the file already exists, we need to remove it!
if [ -f $3/conditional_snp.list ]
then
	rm $3/conditional_snp.list
fi

for chr in {1..22} X
do

	fgrep -f $2 $1/chr${chr}.bimbam | tr "," " " >> $3/conditional_snp.list

done	

transpose.sh $3/conditional_snp.list > $3/conditional_snp_transp.list
