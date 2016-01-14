#!/usr/bin/env bash

IFS=$'\n'
set $(cat metabolic_FVG_parameter_file.txt)
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

pheno=${6} # <- args[[1]]
trait=${4} # <- args[[2]]
covariates=${8} # <- args[[3]]
kinship=${10} # <- args[[4]]
geno=${12} # <- args[[5]]
cohort=${2} # <- args[[6]]
imp_path=${16} # <- args[[7]]
out_path=${18} # <- args[[8]]


for chr in $(seq 20 22)
do
# RUN GWAS analyses using the GWA function
qsub -N "${cohort}_chr${chr}_${trait}" -o "${out_path}/${cohort}_chr${chr}_${trait}.o" \
-e "${out_path}/${cohort}_chr${chr}_${trait}.e" \
-l h_rt=200:00:00 -l vf=10G -wd ${out_path} -q all.q R CMD BATCH '--args '${pheno}' '${trait}' '${covariates}' '${kinship}' '${geno}' '${cohort}' '${imp_path}' '${out_path}'' GWAS_1KG_imputed.R

done

# RUN FORMATTING AND PLOTTING STEPS
bsub -o formatting.log -J "formatting" -w "ended(P_$2_$4*)" \
-M6000000 -R"select[mem>6000] rusage[mem=6000]" -q long \
R CMD BATCH '--args '$2' '$4'' formatting.R

