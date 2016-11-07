#!/usr/bin/env bash
set -e

# cohort="CARL"
# #this can be a different file for logistic and linear analyses
# trait_list="/home/cocca/analyses/CKDgen/ALL_traits_LINEAR.txt"
# pheno_file="/home/cocca/analyses/CKDgen/CARL/pheno/out-CARL/CARL_final_R4.ped"
# imp_path="/netapp02/data/imputation/HRC/${cohort}"
# out_path="/home/cocca/analyses/CKDgen/${cohort}/"

cohort=$1 # <- args[[6]]
trait_list=$2 # <- args[[2]]
pheno_file=$3 # <- args[[1]]
imp_path=$4 # <- args[[7]]
out_path=$5 # <- args[[8]]
mode=$6 # <- args[[8]]

source ~/scripts/bash_scripts/SGE_script_create_function

DIR=/netapp/nfs/softwares/EPACTS/bin
VCF=${imp_path}
PED=${pheno_file}
kinfolder=${out_path}/kinship/${cohort}
mkdir --p ${out_path}


case $mode in
	KINSHIP ) 
		mkdir --p ${out_path}/kinship
		echo "Creating Empirical Kinship matrix..."
		#here we just need to specify the first vcf file
		echo "${DIR}/epacts make-kin --vcf ${VCF}/chr1.dose.vcf.gz --min-maf 0.01 -sepchr --out ${kinfolder}.single.q.emmax.kinf --run 5" | qsub -N KIN_${cohort} -o \$JOB_ID_KIN_${cohort}.o -e \$JOB_ID_KIN_${cohort}.e -V -l h_vmem=2G -cwd
	;;
	LINEAR )

		while read trait
		do

		for chr in $(seq 1 22)
		do
		# RUN GWAS analyses using the GWA function
		# we need to create the script, than we'll submit it
		OUT=${out_path}/${mode}/${trait}/${chr}
		mkdir -p ${out_path}/${mode}/${trait}
		# sge_script_create "${cohort}_chr${chr}_${trait}" "${out_path}/${cohort}_chr${chr}_${trait}.o" "${out_path}/${cohort}_chr${chr}_${trait}.e" ${out_path} R CMD BATCH \'--args ${pheno} ${trait} ${covariates} ${kinship} ${geno} ${cohort} ${chr} ${imp_path}\' ~/scripts/r_scripts/GWAS_1KG_imputed.R ${out_path}/MetS_score_analysis_chr${chr}.Rout > ${out_path}/MetS_score_analysis_chr${chr}.sh
		# here we need to work by chromosome
		sge_script_create "${cohort}_chr${chr}_${trait}" "${out_path}/${cohort}_chr${chr}_${trait}.o" "${out_path}/${cohort}_chr${chr}_${trait}.e" ${out_path}/${trait}  ${DIR}/epacts single --vcf ${VCF}/chr${chr}.dose.vcf.gz --ped ${PED} --chr ${chr} --pheno ${trait} --unit 10000 --test q.emmax --out ${OUT}.single.q.emmax --kinf ${kinfolder}.single.q.emmax.kinf --run 1 > ${out_path}/${mode}/${trait}/CKDGEN_R4_${trait}_chr${chr}.sh

		chmod ug+x ${out_path}/${mode}/${trait}/CKDGEN_R4_${trait}_chr${chr}.sh

		# echo "${out_path}/${mode}/${trait}/CKDGEN_R4_${trait}_chr${chr}.sh" | qsub -N "${cohort}_chr${chr}_${trait}" -o "${out_path}/${cohort}_chr${chr}_${trait}.o" -e "${out_path}/${cohort}_chr${chr}_${trait}.e" -cwd -V -hold_jid KIN_${cohort} -l h_vmem=3G
		echo "${out_path}/${mode}/${trait}/CKDGEN_R4_${trait}_chr${chr}.sh" | qsub -N "${cohort}_chr${chr}_${trait}" -o "${out_path}/${cohort}_chr${chr}_${trait}.o" -e "${out_path}/${cohort}_chr${chr}_${trait}.e" -cwd -V -l h_vmem=3G
		done
		done < <(cat ${trait_list})

	;;
	LOGISTIC )
		btest=$7
		while read trait
		do

		for chr in $(seq 1 22)
		do
		# RUN GWAS analyses using the GWA function
		# we need to create the script, than we'll submit it
		OUT=${out_path}/${mode}/${trait}/${chr}
		mkdir -p ${out_path}/${mode}/${trait}
		# here we need to work by chromosome
		# echo "Running LOGISTIC regression model test..."
		case $trait in
			Gout_men )
				sge_script_create "${cohort}_chr${chr}_${trait}" "${out_path}/${cohort}_chr${chr}_${trait}.o" "${out_path}/${cohort}_chr${chr}_${trait}.e" ${out_path}/${trait}  ${DIR}/epacts single --vcf ${VCF}/chr${chr}.dose.vcf.gz --ped ${PED} --chr ${chr} --pheno ${trait} --cov AGE --cov PC1 --cov PC2 --cov PC3 --cov PC4 --cov PC5 --cov PC6 --cov PC7 --cov PC8 --cov PC9 --cov PC10 --unit 10000 --test ${btest} --out ${OUT}.single.b.wald  --run 1 > ${out_path}/${mode}/${trait}/CKDGEN_R4_${trait}_chr${chr}.sh
				;;
			Gout_women )
				sge_script_create "${cohort}_chr${chr}_${trait}" "${out_path}/${cohort}_chr${chr}_${trait}.o" "${out_path}/${cohort}_chr${chr}_${trait}.e" ${out_path}/${trait}  ${DIR}/epacts single --vcf ${VCF}/chr${chr}.dose.vcf.gz --ped ${PED} --chr ${chr} --pheno ${trait} --cov AGE --cov PC1 --cov PC2 --cov PC3 --cov PC4 --cov PC5 --cov PC6 --cov PC7 --cov PC8 --cov PC9 --cov PC10 --unit 10000 --test ${btest} --out ${OUT}.single.b.wald  --run 1 > ${out_path}/${mode}/${trait}/CKDGEN_R4_${trait}_chr${chr}.sh
				;;
				* )
				sge_script_create "${cohort}_chr${chr}_${trait}" "${out_path}/${cohort}_chr${chr}_${trait}.o" "${out_path}/${cohort}_chr${chr}_${trait}.e" ${out_path}/${trait}  ${DIR}/epacts single --vcf ${VCF}/chr${chr}.dose.vcf.gz --ped ${PED} --chr ${chr} --pheno ${trait} --cov AGE --cov SEX --cov PC1 --cov PC2 --cov PC3 --cov PC4 --cov PC5 --cov PC6 --cov PC7 --cov PC8 --cov PC9 --cov PC10 --unit 10000 --test ${btest} --out ${OUT}.single.b.wald  --run 1 > ${out_path}/${mode}/${trait}/CKDGEN_R4_${trait}_chr${chr}.sh
				;;
		esac
		chmod ug+x ${out_path}/${mode}/${trait}/CKDGEN_R4_${trait}_chr${chr}.sh

		echo "${out_path}/${mode}/${trait}/CKDGEN_R4_${trait}_chr${chr}.sh" | qsub -N "${cohort}_chr${chr}_${trait}" -o "${out_path}/${cohort}_chr${chr}_${trait}.o" -e "${out_path}/${cohort}_chr${chr}_${trait}.e" -cwd -V -l h_vmem=3G
		done
		done < <(cat ${trait_list})
	;;
esac



