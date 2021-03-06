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
			gene_name=`echo ${line}| cut -f 1 -d " "`
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

echo "$FORMAT,$TYPE"
echo "Analyzing region => ${chr}:${start}-${end} ,"
echo "Length: ${gene_length}."
echo "Var number: ${var_num}."
echo "Gene: ${gene_name},"
echo "NExons: ${exon_count}."

#we don't want to keep all regions with 0 variants found, so...
if [[ $var_num -gt 0 ]]; then
#first a count of HET sites
n_het=`egrep "^PSC" ${OUT_F}/${gene_name}/WES.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.stats | awk '$6!=0'| awk 'END{print NR}'`

#than a count of altHOM sites
n_althom=`egrep "^PSC" ${OUT_F}/${gene_name}/WES.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.stats | awk '$5!=0'| awk 'END{print NR}'`

#add the count of samples with either het or hom mutations
n_tot_mut=`egrep "^PSC" ${OUT_F}/${gene_name}/WES.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.stats | awk '{print $5+$6}'|awk '$1!=0'| awk 'END{print NR}'`

#check indels number
n_indels=`egrep "^PSC" ${OUT_F}/${gene_name}/WES.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.stats | awk '$9!=0'| awk 'END{print NR}'`

#mean coverage in region
mean_cov=`egrep "^PSC" ${OUT_F}/${gene_name}/WES.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.stats| awk '{sum += $10} END {print sum / NR }'`

#this is the only number that makes sense, because we're divinding the number of samples with het or hom mutation, by the total number of samples
perc_het_samples=`printf "%f\n" $( bc -l <<< "${n_het}/${all_inds}")`
perc_althom_samples=`printf "%f\n" $( bc -l <<< "${n_althom}/${all_inds}")`
perc_tot_mut_samples=`printf "%f\n" $( bc -l <<< "${n_tot_mut}/${all_inds}")`
perc_indels_samples=`printf "%f\n" $( bc -l <<< "${n_indels}/${all_inds}")`

#we should also normalize by the number of variants in our gene
#but this should be done on the sample frequency
perc_het_nsites=`printf "%f\n" $( bc -l <<< "${n_het}/${var_num}")`
perc_althom_nsites=`printf "%f\n" $( bc -l <<< "${n_althom}/${var_num}")`
perc_het_length=`printf "%f\n" $( bc -l <<< "${n_het}/${gene_length}")`
perc_althom_length=`printf "%f\n" $( bc -l <<< "${n_althom}/${gene_length}")`


case ${TYPE} in
	SNP )
	#now print a resume line for this gene
	(echo "ID CHR START END EXON_COUNT GENE_LENGTH VARIANT_NUMBER TOT_SAMPLES HET_SAMPLES ALT_HOM_SAMPLES TOT_MUT_SAMPLES FREQ_HET_BY_SAMPLE FREQ_ALT_HOM_BY_SAMPLE FREQ_TOT_MUT_BY_SAMPLE FREQ_HET_BY_SITE FREQ_ALT_HOM_BY_SITE FREQ_HET_BY_LENGTH FREQ_ALT_HOM_BY_LENGTH GENE MEAN_COVERAGE"
	echo ${gene_name} ${chr} ${start} ${end} ${exon_count} ${gene_length} ${var_num} ${all_inds} ${n_het} ${n_althom} ${n_tot_mut} ${perc_het_samples} ${perc_althom_samples} ${perc_tot_mut_samples} ${perc_het_nsites} ${perc_althom_nsites} ${perc_het_length} ${perc_althom_length} ${gene} ${mean_cov}) > ${OUT_F}/${gene_name}/WES.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.stats.resume
		;;
	INDEL )
	#now print a resume line for this gene
	(echo "ID CHR START END EXON_COUNT GENE_LENGTH VARIANT_NUMBER TOT_SAMPLES HET_SAMPLES ALT_HOM_SAMPLES FREQ_HET_BY_SAMPLE FREQ_ALT_HOM_BY_SAMPLE FREQ_HET_BY_SITE FREQ_ALT_HOM_BY_SITE FREQ_HET_BY_LENGTH FREQ_ALT_HOM_BY_LENGTH GENE MEAN_COVERAGE"
	echo ${gene_name} ${chr} ${start} ${end} ${exon_count} ${gene_length} ${var_num} ${all_inds} ${n_indels} ${n_althom} ${perc_indels_samples} ${perc_althom_samples} ${perc_het_nsites} ${perc_althom_nsites} ${perc_het_length} ${perc_althom_length} ${gene} ${mean_cov}) > ${OUT_F}/${gene_name}/WES.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.stats.resume

		;;
esac
else
	echo "Empty region!!Removing useless files!!"
	rm -rf ${OUT_F}/${gene_name}/
fi

done < <(zcat $REG_F)
