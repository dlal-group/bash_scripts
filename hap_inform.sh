#!/usr/bin/env bash
#20/01/2017
# Pipeline to extract most informative haplotypes from 1000G phase 3 data

TGP_input=$1
base_out=$2
maf=$3
sd=$4
subset_sample=$5


# TGP_input="/netapp/nfs/resources/1000GP_phase3/vcf/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
# base_out="/home/cocca/analyses/michelangelo"
# maf=0.5
# sd=0.05
# subset_sample="/home/cocca/analyses/michelangelo/EUR_samples_phase3.txt"
#add subset Asian and whole TGP3

basename_out=`basename ${TGP_input}`

ad_maf1=`echo "$maf $sd" | awk '{print $1-$2}'`
# ad_maf2=$[maf + sd]
#extract MAF filtered data
bcftools view -i"MAF>=${ad_maf1}" -v snps -S ${subset_sample} ${TGP_input} -O z -o ${base_out}/${basename_out}.${maf}.vcf.gz 

#Extract snp list based on distance: we want to maximize the number of sites in a 100bp window