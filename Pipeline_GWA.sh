#!/usr/bin/env bash

IFS=$'\n'
set $(cat ${1})
# set $(cat parameter_file.txt)

# Cohort's name $2 = FVG
# TRAIT $4 = MetS_score
# Pheno file GenABEL's input $6 = /home/cocca/analyses/MetabolicSyndrome/FVG/fvg_all_metabolic_ALL_MetS_score.csv
# Covariates $8 = age
# Genomic kinship $10 = /nfs/servizio/FVG.kinship
# Genotypes file $12 = /home/cocca/analyses/MetabolicSyndrome/FVG/FVG_out
# Map file $14 = /nfs/1000G/FVG/dose/FVG_1000G_
# Info file $16 = /nfs/1000G/FVG/prob/FVG_1000G_
# Dose file $18 = /nfs/1000G/FVG/dose/FVG_1000G_
# Output folder $20 

# pheno="/home/cocca/analyses/MetabolicSyndrome/FVG/fvg_all_metabolic_ALL_MetS_score.csv" # <- args[[1]]
# trait="MetS_score" # <- args[[2]]
# covariates="age" # <- args[[3]]
# kinship="/nfs/servizio/FVG.kinship" # <- args[[4]]
# geno="/home/cocca/analyses/MetabolicSyndrome/FVG/FVG_out" # <- args[[5]]
# cohort="FVG" # <- args[[6]]
# imp_path="/nfs/1000G/FVG/prob/FVG_1000G" # <- args[[7]]
# out_path="/home/cocca/analyses/MetabolicSyndrome/FVG/" # <- args[[8]]


pheno=${6} # <- args[[1]]
trait=${4} # <- args[[2]]
covariates=${8} # <- args[[3]]
kinship=${10} # <- args[[4]]
geno=${12} # <- args[[5]]
cohort=${2} # <- args[[6]]
imp_path=${16} # <- args[[7]]
out_path=${20} # <- args[[8]]

echo -e "pheno = ${pheno}\n
trait = ${trait}\n
covariates = ${covariates}\n
kinship = ${kinship}\n
geno = ${geno}\n
cohort = ${cohort}\n
imp_path = ${imp_path}\n
out_path = ${out_path}\n
"
source ~/scripts/bash_scripts/SGE_script_create_function

for chr in $(seq 1 22)
do
# RUN GWAS analyses using the GWA function
# we need to create the script, than we'll submit it
mkdir -p ${out_path}
# sge_script_create "${cohort}_chr${chr}_${trait}" "${out_path}/${cohort}_chr${chr}_${trait}.o" "${out_path}/${cohort}_chr${chr}_${trait}.e" ${out_path} R CMD BATCH \'--args ${pheno} ${trait} ${covariates} ${kinship} ${geno} ${cohort} ${chr} ${imp_path}\' ~/scripts/r_scripts/GWAS_1KG_imputed.R ${out_path}/MetS_score_analysis_chr${chr}.Rout > ${out_path}/MetS_score_analysis_chr${chr}.sh
sge_script_create "${cohort}_chr${chr}_${trait}" "${out_path}/${cohort}_chr${chr}_${trait}.o" "${out_path}/${cohort}_chr${chr}_${trait}.e" ${out_path} R CMD BATCH \'--args ${pheno} ${trait} ${covariates} ${kinship} ${geno} ${cohort} ${chr} ${imp_path}\' ~/scripts/r_scripts/GWAS_1000G_phd.R ${out_path}/1000G_PHD_${trait}_chr${chr}.Rout > ${out_path}/1000G_PHD_${trait}_chr${chr}.sh
chmod ug+x ${out_path}/1000G_PHD_${trait}_chr${chr}.sh

qsub ${out_path}/1000G_PHD_${trait}_chr${chr}.sh
done

