#!/usr/bin/env bash

# This is the runner file run by SGE
# Arguments: runner.sh filelist
# Environment variables: SGE_TASK_ID
# a_size=`wc -l chr${chr}_command.list| cut -f 1 -d " "`;echo "~/scripts/bash_scripts/ja_runner_par_TRST.sh -l $imputedir/chr${chr}_command.list"|qsub -t 1-${a_size} -o ${imputedir}/chr${chr}_\$JOB_ID_\$TASK_ID.log -e ${imputedir}/chr${chr}_\$JOB_ID_\$TASK_ID.e -V -N ${pop}_chr${chr} -l h_vmem=${m}
set -e

file=`sed -n "${SGE_TASK_ID}p" $1`

#11/09/2016
#remove duplicate lines from vcf files
base_dir=`dirname ${file}`
file_name=`basename ${file}`
mkdir -p ${base_dir}/11092016_ANN
mkdir -p ${base_dir}/11092016_ANN/TAB_SNP
mkdir -p ${base_dir}/11092016_ANN/TAB_INDEL

# echo "Cleaning possible duplicate rows in annotated vcfs...."

# (bcftools view -h ${file};bcftools view -H ${file}| uniq)|bgzip -c > ${base_dir}/11092016_ANN/${file_name}
# tabix -f -p vcf ${base_dir}/11092016_ANN/${file_name}

# #remove also eventually duplicate lines from CADD annotation files
# echo "Cleaning possible duplicate rows in annotation files...."
# zcat ${base_dir}/TAB/${file_name}.scores.tsv.gz| uniq | gzip -c > ${base_dir}/11092016_ANN/TAB_SNP/${file_name}.scores.tsv.gz
# zcat ${base_dir}/TAB_INDEL/${file_name}.scores.tsv.gz| uniq | gzip -c > ${base_dir}/11092016_ANN/TAB_INDEL/${file_name}.scores.tsv.gz

#Preparing CADD tables to annotate vcf files
#SNP section

#we need a table to correctly sort multiple alleles
#need to have CHR POS REF ALT, but only for multiallelic sites
chr=`echo ${file_name%%.*}`

#here we work for indels and snp
for type in snp indel
do

U_TYPE=`echo ${type^^}`

echo "working on ${U_TYPE}..."

bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" -i"TYPE='${type}'" ${base_dir}/11092016_ANN/${file_name} | awk '{if($4~",") print $0}' > ${base_dir}/11092016_ANN/TAB_${U_TYPE}/${chr}.MULTI_${U_TYPE}.list

echo "chr ${chr} .."
echo "Preparing annotation file...."
case ${type} in
	snp )
		zcat ${base_dir}/11092016_ANN/TAB_${U_TYPE}/${file_name}.scores.tsv.gz| tail -n+2 | awk '{print $1,$2,$3,$5,$(NF-1),$NF}' | prepare_annots.py ${base_dir}/11092016_ANN/TAB_${U_TYPE}/${chr}.MULTI_${U_TYPE}.list | sort -g -k2,2 | bgzip -c > ${base_dir}/11092016_ANN/TAB_${U_TYPE}/${chr}.formatted.CADD.tab.gz
	;;
	indel)
		zcat ${base_dir}/11092016_ANN/TAB_${U_TYPE}/${file_name}.scores.tsv.gz| tail -n+2 | awk '{print $1,$2,$3,$4,$(NF-1),$NF}' | prepare_annots.py ${base_dir}/11092016_ANN/TAB_${U_TYPE}/${chr}.MULTI_${U_TYPE}.list | sort -g -k2,2 | bgzip -c > ${base_dir}/11092016_ANN/TAB_${U_TYPE}/${chr}.formatted.CADD.tab.gz
	;;	
esac

tabix -f -s 1 -b 2 -e 2 ${base_dir}/11092016_ANN/TAB_${U_TYPE}/${chr}.formatted.CADD.tab.gz

mkdir -p ${base_dir}/11092016_ANN/11092016_CADD_ANNOT
echo "Annotate ${U_TYPE} in vcf...."
bcftools view -v ${type}s  ${base_dir}/11092016_ANN/${file_name} | bcftools annotate -a ${base_dir}/11092016_ANN/TAB_${U_TYPE}/${chr}.formatted.CADD.tab.gz -c CHROM,POS,REF,ALT,CADD_RAW,CADD_PHRED -h /netapp/dati/INGI_WGS/18112015/CADD_header.txt -O z -o ${base_dir}/11092016_ANN/11092016_CADD_ANNOT/${chr}.${U_TYPE}.vcf.gz
tabix -f -p vcf ${base_dir}/11092016_ANN/11092016_CADD_ANNOT/${chr}.${U_TYPE}.vcf.gz
echo "Chromosome ${chr} done."

done

#now we join indels and snps back together
echo "Join INDEL and SNP file for chr${chr} ..."
bcftools concat ${base_dir}/11092016_ANN/11092016_CADD_ANNOT/${chr}.SNP.vcf.gz ${base_dir}/11092016_ANN/11092016_CADD_ANNOT/${chr}.INDEL.vcf.gz -O z -o ${base_dir}/11092016_ANN/11092016_CADD_ANNOT/${chr}.JOINT.vcf.gz

echo "Sort by position the JOINT file for chr${chr} ..."
(bcftools view -h ${base_dir}/11092016_ANN/11092016_CADD_ANNOT/${chr}.JOINT.vcf.gz;bcftools view -H ${base_dir}/11092016_ANN/11092016_CADD_ANNOT/${chr}.JOINT.vcf.gz | sort -g -k2,2 ) | bgzip -c > ${base_dir}/11092016_ANN/11092016_CADD_ANNOT/${chr}.vcf.gz
tabix -p vcf ${base_dir}/11092016_ANN/11092016_CADD_ANNOT/${chr}.vcf.gz
