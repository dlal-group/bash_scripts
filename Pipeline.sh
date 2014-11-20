#!/usr/local/bin/bash

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

IFS=$'\n'
set $(cat parameter_file.txt)

# RUN GenABEL
bsub -J "G_$2_$4" -o GenABEL.log  -M6000000 -R"select[mem>6000] rusage[mem=6000]" \
-q long R CMD BATCH '--args '$6' '$4' '$8' '${10}' '${12}'' GenABEL.R

# RUN ProbABEL
for i in $(seq 1 22)
do
bsub -J "P_$2_$4_chr${i}" -w "ended(G_$2_$4)" -o chr${i}.log -M7000000 -R"select[mem>7000] \
rusage[mem=7000]" -q long /nfs/users/nfs_g/gp5/ProbABEL/bin/palinear --mmscore varcovar.mat --pheno res.pheno \
--chrom ${i} --map ${14}_chr${i}.map \
--info ${16}mach${i}Chr.machinfo \
--dose ${18}/chr${i}/mach${i}chr.fvi --out chr${i}.palinear 
done

# RUN FORMATTING AND PLOTTING STEPS
bsub -o formatting.log -J "formatting" -w "ended(P_$2_$4*)" \
-M6000000 -R"select[mem>6000] rusage[mem=6000]" -q long \
R CMD BATCH '--args '$2' '$4'' formatting.R


