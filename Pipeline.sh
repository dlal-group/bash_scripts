#!/usr/local/bin/bash

IFS=$'\n'
set $(cat metabolic_FVG_parameter_file.txt)
# set $(cat parameter_file.txt)

# Cohort's name $2 = FVG
# TRAIT $4 = MetS_score
# Pheno file GenABEL's input $6 = /home/cocca/analyses/MetabolicSyndrome/FVG/fvg_all_metabolic_ALL_MetS_score.csv
#Covariates $8 = age
# Genomic kinship $10 = /nfs/servizio/FVG.kinship
# Genotypes file $12 = /home/cocca/analyses/MetabolicSyndrome/FVG/FVG_out
# Map file $14 = /nfs/1000G/FVG/dose/FVG_1000G_
# Info file $16 = /nfs/1000G/FVG/dose/FVG_1000G_
# Dose file $18 = /nfs/1000G/FVG/dose/FVG_1000G_


# INSERT COHORT'S NAME
#echo -e "Cohort's name: \c "
#read word1
#echo $word1

#while read line
#do
#echo -e $line
#$VAR=$line
#echo -e $VAR
#done <param.txt



# RUN GWAS analyses using the GWA function
bsub -J "G_$2_$4" -o GenABEL.log  -M6000000 -R"select[mem>6000] rusage[mem=6000]" \
-q long R CMD BATCH '--args '${6}' '${4}' '${8}' '${10}' '${12}' '${2}'' GWAS_1KG_imputed.R

# RUN ProbABEL
for i in $(seq 1 22)
do
bsub -J "P_$2_$4_chr${i}" -w "ended(G_$2_$4)" -o chr${i}.log -M7000000 -R"select[mem>7000] \
rusage[mem=7000]" -q long /home/cocca/softwares/probabel/bin/palinear --mmscore varcovar.mat --pheno res.pheno \
--chrom ${i} --map ${14}_chr${i}.map \
--info ${16}_chr${i}.info \
--dose ${18}_chr${i}.fvi --out chr${i}.palinear 
done

# RUN FORMATTING AND PLOTTING STEPS
bsub -o formatting.log -J "formatting" -w "ended(P_$2_$4*)" \
-M6000000 -R"select[mem>6000] rusage[mem=6000]" -q long \
R CMD BATCH '--args '$2' '$4'' formatting.R


/home/cocca/softwares/probabel/bin/palinear --mmscore varcovar.mat --pheno res.pheno --chrom 22 --map /nfs/1000G/FVG/dose/FVG_1000G_chr22.map --info /nfs/1000G/FVG/prob/FVG_1000G_chr22.info --dose /nfs/1000G/FVG/dose/FVG_1000G_chr22.fvi --out chr22.palinear 