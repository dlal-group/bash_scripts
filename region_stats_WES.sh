#!/bin/bash

#script to extract stats for all samples for a single region from annotated vcf file
IN_VCF=$1 #input vcf folder
REG_F=$2 #region file
OUT_F=$3 #output folder
FORMAT=$4 #based on format of our file
TYPE=$5 #variation type : SNP / INDEL

while read line
do
	case ${FORMAT} in
		UCSC )
			chr=`echo ${line}| cut -f 2 -d " "| sed 's/chr//g'`
			start=`echo ${line}| cut -f 4 -d " "`
			end=`echo ${line}| cut -f 5 -d " "`
			gene=`echo ${line}| cut -f 1 -d " "`
			exon_count=`echo ${line}| cut -f 8 -d " "`
			;;
		GENCODE )
			chr=`echo ${line}| cut -f 3 -d " "|sed 's/chr//g'`
			start=`echo ${line}| cut -f 5 -d " "`
			end=`echo ${line}| cut -f 6 -d " "`
			gene_name=`echo ${line}| cut -f 2 -d " "`
			gene=`echo ${line}| cut -f 2 -d " "`
			exon_count=`echo ${line}| cut -f 9 -d " "`
			;;
		GENCODE19 )
			chr=`echo ${line}| cut -f 1 -d " "|sed 's/chr//g'`
			start=`echo ${line}| cut -f 4 -d " "`
			end=`echo ${line}| cut -f 5 -d " "`
			gene_name=`echo ${line}| cut -f 9 | cut -f 5 -d ";" | awk '{print $2}'|sed 's/^"\(.*\)"$/\1/'`
			gene=`echo ${line}| cut -f 9 | cut -f 5 -d ";" | awk '{print $2}'|sed 's/^"\(.*\)"$/\1/'`
			exon_count="NA"
			;;
		TRANSCRIPT19 )
			chr=`echo ${line}| cut -f 1 -d " "|sed 's/chr//g'`
			start=`echo ${line}| cut -f 4 -d " "`
			end=`echo ${line}| cut -f 5 -d " "`
			gene_name=`echo ${line}| cut -f 9 | cut -f 2 -d ";" | awk '{print $2}'|sed 's/^"\(.*\)"$/\1/'`
			gene=`echo ${line}| cut -f 9 | cut -f 5 -d ";" | awk '{print $2}'|sed 's/^"\(.*\)"$/\1/'`
			exon_count="NA"
			;;
	esac
			gene_length=$[end - start]

echo "$FORMAT,$TYPE"
echo "Analyzing region => ${chr}:${start}-${end} ,"
echo "Length: ${gene_length}."
echo "Gene: ${gene_name},"
echo "NExons: ${exon_count}."

mkdir -p ${OUT_F}/${gene_name}
OUT_VCF=${OUT_F}/${gene_name}/All.multisampleinitial.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.recalibrated.filtered.vcf.gz


# we need to read from our region file and save the extracted region, than we're going to extract stats for that region
bcftools view ${IN_VCF} -r ${chr}:${start}-${end} -O z -o ${OUT_VCF}
tabix -f -p vcf ${OUT_VCF}
bcftools stats -s - ${OUT_VCF} > ${OUT_F}/${gene_name}/WES.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.stats
			
#we need to extract the Per Sample Count section of our stat
#and than calculate how many samples have nNonRefHom mutations and how many nHets
#and how many singletons and add also gene informations (size,region, number of sites)

var_num=`egrep "^SN" ${OUT_F}/${gene_name}/WES.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.stats| fgrep "number of records"| cut -f 4`
all_inds=`egrep "^SN" ${OUT_F}/${gene_name}/WES.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.stats| fgrep "number of samples"| cut -f 4`
#we don't want to keep all regions with 0 variants found, so...
if [[ $var_num -gt 0 ]]; then
#first a count of HET sites
n_het=`egrep "^PSC" ${OUT_F}/${gene_name}/WES.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.stats | awk '$6!=0'| awk 'END{print NR}'`

#than a count of altHOM sites
n_althom=`egrep "^PSC" ${OUT_F}/${gene_name}/WES.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.stats | awk '$5!=0'| awk 'END{print NR}'`

#this is the only number that makes sense, because we're divinding the number of samples with het or hom mutation, by the total number of samples
perc_het_samples=$( bc -l <<< "${n_het}/${all_inds}")
perc_althom_samples=$( bc -l <<< "${n_althom}/${all_inds}")

#we should also normalize by the number of variants in our gene
#but this should be done on the sample frequency
perc_het_nsites=$( bc -l <<< "${n_het}/${var_num}")
perc_althom_nsites=$( bc -l <<< "${n_althom}/${var_num}")
perc_het_length=$( bc -l <<< "${n_het}/${gene_length}")
perc_althom_length=$( bc -l <<< "${n_althom}/${gene_length}")


#now print a resume line for this gene
(echo "ID CHR START END EXON_COUNT GENE_LENGTH VARIANT_NUMBER TOT_SAMPLES HET_SAMPLES ALT_HOM_SAMPLES FREQ_HET_BY_SAMPLE FREQ_ALT_HOM_BY_SAMPLE FREQ_HET_BY_SITE FREQ_ALT_HOM_BY_SITE FREQ_HET_BY_LENGTH FREQ_ALT_HOM_BY_LENGTH GENE"
echo ${gene_name} ${chr} ${start} ${end} ${exon_count} ${gene_length} ${var_num} ${all_inds} ${n_het} ${n_althom} ${perc_het_samples} ${perc_althom_samples} ${perc_het_nsites} ${perc_althom_nsites} ${perc_het_length} ${perc_althom_length} ${gene}) > ${OUT_F}/${gene_name}/WES.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.stats.resume
else
	echo "Empty region!!Removing useless files!!"
	rm -rf ${OUT_F}/${gene_name}/
fi

done < <(zcat $REG_F)
