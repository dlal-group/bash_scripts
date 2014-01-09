#!/usr/local/bin/bash

# This is the runner file run by bsub
# Arguments: runner.sh filelist
# Environment variables: LSB_JOBINDEX

file=`sed -n "${LSB_JOBINDEX}p" $2`

# script=$1
# uncomment this if you need to work with chr as jobindex
# if [[ ${LSB_JOBINDEX} -eq 23 ]]
# then
# 	chr="X"
# else
# 	chr=${LSB_JOBINDEX}
# fi
# chr=${LSB_JOBINDEX}
# chr=$2
# cohort=$1
# trait=$1
outpath=$1
# bash $script $file $3 $4 $5
# /software/vertres/codebase/scripts/bamcheck -c 1,50,1 -d /lustre/scratch113/projects/fvg_seq/F12HPCEUHK0358_HUMpngR/BRIDGED_BAMS/${file}
# zcat ${file} | sed 's/,rs/|rs/g' | awk -v chr=${LSB_JOBINDEX} '{ snp=(NF-5)/3; if($2 ~/^rs/) s=$2;else s="NA"; printf "chr"chr":"$3"-"s"-"$4"-"$5"," $4 "," $5; for(i=1; i<=snp; i++) printf "," $(i*3+3)*2+$(i*3+4); printf "\n" }' > chr${LSB_JOBINDEX}.bimbam
#Without and with permutation
# plink --noweb --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bed --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bim.COPY --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.fam --allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/${cohort}_samples.keeplist --pheno /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/${cohort}_samples_pheno.centre --assoc --out chr${chr}.txt
#plink --noweb --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bed --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bim.COPY --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.fam --allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/${cohort}_samples.keeplist --pheno /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/${cohort}_samples_pheno.centre --assoc --perm --out chr${chr}.txt
# plink --noweb --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bed --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bim.COPY --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.fam --allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/${cohort}_samples.keeplist --pheno /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/${cohort}_samples_pheno.centre --assoc --adjust --out chr${chr}.txt
# plink --noweb --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bed --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bim.COPY --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.fam --allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/3621_${cohort}_samples.keeplist --pheno /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/${cohort}_samples_pheno.centre --assoc --adjust --out chr${chr}.txt
#plink --noweb --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bed --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bim.COPY --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.fam --allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/${cohort}_samples.keeplist --pheno /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/${cohort}_samples_pheno.centre --maf 0.05 --assoc --adjust --out chr${chr}.txt
#sed 's/ \+/ /g' chr${chr}.txt.assoc | sed 's/^ //g' | tr " " "\t" | gzip -c > chr${chr}.txt.assoc.gz
#plink --noweb --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bed --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bim.COPY --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.fam --allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/${cohort}_samples.keeplist --pheno /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/${cohort}_samples_pheno.centre --max-maf 0.05 --assoc --adjust --out chr${chr}.txt
# plink --noweb --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bed --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bim.COPY --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.fam --allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/${cohort}_samples.keeplist --freq --out chr${chr}.txt
# plink --noweb --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bed --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bim --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.fam --allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/${cohort}_samples.keeplist --freq --out chr${chr}.pruned.txt
#unfiltered set
#plink --noweb --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bed --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bim.COPY --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.fam --allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/3621_${cohort}_samples.keeplist --freq --out chr${chr}.freq
#mkdir -p REMOVELIST
# cat ${file}| sed 's/ \+/ /g'|sed 's/^ //g' | sed 's/ $//g' | sed 's/ /\t/g' |awk '$9<5e-8' | cut -f 2 > REMOVELIST/${file}.$1.5e-8.removelist
# cat ${file}| sed 's/ \+/ /g'|sed 's/^ //g' | sed 's/ $//g' | sed 's/ /\t/g' |awk '$9<5e-6' | cut -f 2 > REMOVELIST/${file}.$1.5e-6.removelist
# cat ${file}| sed 's/ \+/ /g'|sed 's/^ //g' | sed 's/ $//g' | sed 's/ /\t/g' |awk '$9<1e-5' | cut -f 2 > REMOVELIST/${file}.$1.1e-5.removelist
#cat ${file}| sed 's/ \+/ /g'|sed 's/^ //g' | sed 's/ $//g' | sed 's/ /\t/g' |awk '$9<1e-2' | cut -f 2 > REMOVELIST/${file}.$1.1e-2.removelist
# fgrep -v -w -f ~/UK10K/users/jh21/imputed/fvg/uk10k1kg.shapeit/master_stats/rare_variants.list /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/GLYCEMIC/output/${trait}.chr${chr}.tab.assoc.txt > ${trait}.chr${chr}.gt001.tab.assoc.txt

#CC stratified by population - maf > 1%
#plink --noweb \
#--bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bed \
#--bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bim.COPY \
#--fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.fam \
#--allow-no-sex \
#--keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples.keeplist \
#--pheno /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples_pheno.centre \
#--maf 0.01 \
#--logistic \
##--assoc \
##--mh \
#--covar /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples_pheno.pop \
##--within /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples_pheno.pop \
#--adjust --out chr${chr}.txt
#calculate freq for the single chr
#mkdir -p /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/UNFILTERED/CHR${chr}

#plink --noweb \
#--bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
#--allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples.keeplist \
#--chr ${chr} \
#--freq \
#--out /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/UNFILTERED/CHR${chr}/chr${chr}_freq

#CC stratified by population - all sites
# plink --noweb \
# --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bed \
# --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bim.COPY \
# --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.fam \
# --allow-no-sex \
# --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples.keeplist \
# --pheno /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples_pheno.centre \
# --logistic \
# --covar /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples_pheno.pop \
# --adjust --out chr${chr}.txt

#CC logistic alspac and twins - maf >= 1%
# plink --noweb \
# --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bed \
# --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bim.COPY \
# --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.fam \
# --allow-no-sex \
# --maf 0.01 \
# --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples.keeplist \
# --pheno /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples_pheno.centre \
# --logistic \
# --adjust --out chr${chr}.txt

# gzip all results files
# gzip ${file}

# generate stats for bamfiles
# samtools_stats.sh ${file} ${outpath}

#merge genotypes
#chr=${file}
#plink --bfile ~/UK10K/users/jh21/imputed/fvg/fvg_370/shapeit/chr${chr} --bmerge ~/UK10K/users/jh21/imputed/fvg/fvg_omni/shapeit/chr${chr}.bed ~/UK10K/users/jh21/imputed/fvg/fvg_omni/shapeit/chr${chr}.bim ~/UK10K/users/jh21/imputed/fvg/fvg_omni/shapeit/chr${chr}.fam --make-bed --out ${outpath}/chr${chr}_merged

#extract sites by region using antitumoral list
chr=`echo ${file} | cut -f 1 -d " "`
start=`echo ${file} | cut -f 2 -d " "`
end=`echo ${file} | cut -f 3 -d " "`
gene=`echo ${file} | cut -f 4 -d " "`
#plink --noweb --bfile /nfs/users/nfs_m/mc14/Work/SANGER/FVG/ANTI_TUMORAL_DRUGS/merged/chr${chr}_merged --chr ${chr} --from-bp ${start} --to-bp ${end} --make-bed --out ${outpath}/chr${chr}_${gene}

#now calculate also the frequencies in FVG
#plink --noweb --bfile /nfs/users/nfs_m/mc14/Work/SANGER/FVG/ANTI_TUMORAL_DRUGS/merged/chr${chr}_merged --chr ${chr} --from-bp ${start} --to-bp ${end} --freq --out ${outpath}/chr${chr}_${gene}_fvgfrq

#now calculate also the frequencies in FVG for males only
plink --noweb --bfile /nfs/users/nfs_m/mc14/Work/SANGER/FVG/ANTI_TUMORAL_DRUGS/merged/chr${chr}_merged --chr ${chr} --from-bp ${start} --to-bp ${end} --filter-males --freq --out ${outpath}/males/chr${chr}_${gene}_fvgfrq

#now calculate also the frequencies in FVG for females only
plink --noweb --bfile /nfs/users/nfs_m/mc14/Work/SANGER/FVG/ANTI_TUMORAL_DRUGS/merged/chr${chr}_merged --chr ${chr} --from-bp ${start} --to-bp ${end} --filter-females --freq --out ${outpath}/females/chr${chr}_${gene}_fvgfrq

