#!/usr/local/bin/bash

# This is the runner file run by bsub
# Arguments: runner.sh filelist
# Environment variables: LSB_JOBINDEX
# mkdir -p LOGS;size=`wc -l result.list|cut -f 1 -d " "`;bsub -J "p_check[1-${size}]" -o "LOGS/%J_p_check.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner.sh result.list
file=`sed -n "${LSB_JOBINDEX}p" $1`
#added for population control
# pop=$2

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
# outpath=$2
# mkdir -p ${outpath}
# grep_file=$3
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
# --allow-no-sex \file1
# file
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

# generate index for bamfiles
# if [ ! -f ${file}.bai ];
# then
#   samtools index ${file}
# fi

#merge genotypes
#chr=${file}
#plink --bfile ~/UK10K/users/jh21/imputed/fvg/fvg_370/shapeit/chr${chr} --bmerge ~/UK10K/users/jh21/imputed/fvg/fvg_omni/shapeit/chr${chr}.bed ~/UK10K/users/jh21/imputed/fvg/fvg_omni/shapeit/chr${chr}.bim ~/UK10K/users/jh21/imputed/fvg/fvg_omni/shapeit/chr${chr}.fam --make-bed --out ${outpath}/chr${chr}_merged

#extract sites by region using antitumoral list
# echo ${file}
# chr=`echo ${file} | cut -f 1 -d " "`
# start=`echo ${file} | cut -f 2 -d " "`
# end=`echo ${file} | cut -f 3 -d " "`
# gene=`echo ${file} | cut -f 4 -d " "`
# genotype_path=$2
# outpath=$3

# mkdir -p ${outpath}/males
# mkdir -p ${outpath}/females

# #create bed files
# plink --noweb --bfile ${genotype_path} --chr ${chr} --from-bp ${start} --to-bp ${end} --make-bed --out ${outpath}/chr${chr}_${gene}

# #now calculate also the frequencies in FVG
# plink --noweb --bfile ${genotype_path} --chr ${chr} --from-bp ${start} --to-bp ${end} --freq --out ${outpath}/chr${chr}_${gene}_fvgfrq

# #now calculate also the frequencies in FVG for males only
# plink --noweb --bfile ${genotype_path} --chr ${chr} --from-bp ${start} --to-bp ${end} --filter-males --freq --out ${outpath}/males/chr${chr}_${gene}_fvgfrq

# #now calculate also the frequencies in FVG for females only
# plink --noweb --bfile ${genotype_path} --chr ${chr} --from-bp ${start} --to-bp ${end} --filter-females --freq --out ${outpath}/females/chr${chr}_${gene}_fvgfrq

#calculate md5sum
# filename=`basename ${file}`
# md5sum ${file} > ${outpath}/${filename}.md5sum

#15/01/2014
#Extract list of snps for Nicola
# (head -1 ${file};fgrep -w -f ${grep_file} ${file}) > ${outpath}/${file}.nicola

# (zcat ${file}| head -1;zcat ${file}|fgrep -w -f ${grep_file}) > ${outpath}/${file}.nicola

#tabix index vcf file
# tabix -p vcf ${file}.vcf.gz

# #split file in chromosomes
# filename=`basename $2`
# tabix -h $2 ${file} | bgzip -c > $3/${filename}.${file}.vcf.gz
# tabix -p vcf $3/$2.${file}.vcf.gz

#Check sex in ba files
# echo "${file}"
# samtools view ${file} Y| cut -f 1 | sort -u | wc -l

#17/02/2014 subset vcf file
#ARGS:
# $2=vcf file path
# $3= outpath
# village_file=`basename ${file}`
# village=${village_file%_*}


# bcftools2 view -s ${file} $2 -O v | vcf-annotate --fill-ICF | bgzip -c > $3.${village}.vcf.gz
# tabix -f -p vcf $3.${village}.vcf.gz

#extract stats by village
# bcftools2 stats -d 0,5000,1 -s - ${file} > ${file}.vcfchk;plot-vcfstats ${file}.vcfchk -p vcf_check/plots/${file}

# extract table info by village
#(echo "CHROM POS ID REF ALT AN AC TGP_AF AMR_AF ASN_AF AFR_AF EUR_AF IMP2 VQSLOD CULPRIT AF MAF MINOR";bcftools2 query ${file} -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\t%INFO/1kg_AF\t%INFO/1kg_AMR_AF\t%INFO/1kg_ASN_AF\t%INFO/1kg_AFR_AF\t%INFO/1kg_EUR_AF\t%INFO/IMP2\t%INFO/VQSLOD\t%INFO/culprit\n" | awk '{if($5 !~ ",") print $0}' | awk '{if($6 != 0)print $0,$7/$6;else print $0,"NA"}' | awk '{if($NF != "NA") {if($NF > 0.5) print $0,1-$NF,$4;else print $0,$NF,$5}else{print $0,"NA","NA"}}') | tr "\t" " " > ${file}.csv
# (echo "CHROM POS ID REF ALT AN AC VQSLOD CULPRIT AF MAF MINOR";bcftools2 query ${file} -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\t%INFO/VQSLOD\t%INFO/culprit\n" | awk '{if($5 !~ ",") print $0}' | awk '{if($6 != 0)print $0,$7/$6;else print $0,"NA"}' | awk '{if($NF != "NA") {if($NF > 0.5) print $0,1-$NF,$4;else print $0,$NF,$5}else{print $0,"NA","NA"}}') | tr "\t" " " > ${file}.csv

#3/3/2014
#split SNPS and indels
# output1=$2
# output2=$3
# mkdir -p ${output1}
# mkdir -p ${output2}

# base=`basename ${file}`
# zcat ${file} | awk 'length($4) != 1 || length($5) != 1'| cut -f 2- -d " " | awk '{print $1,$0}'| tr "\t" " "| gzip -c > ${output1}/${base}
# zcat ${file} | awk 'length($4) == 1 && length($5) == 1'| cut -f 2- -d " " | awk '{print $1,$0}'| tr "\t" " "| gzip -c > ${output2}/${base}

#3/3/2014
#merge genotype imputation files
# output=$2
# mkdir -p ${output}

# file1=`echo ${file} | cut -f 1 -d " "`
# file2=`echo ${file} | cut -f 2 -d " "`
# base1=`basename ${file1}`
# base2=`basename ${file2}`

# samples1=`echo ${file1%.*}_samples`
# samples2=`echo ${file2%.*}_samples`
# gtool -M --g ${file1} ${file2} --s ${samples1} ${samples2} --log ${output}/${base1}_${base2}.log --og ${output}/${base1}_${base2}.gen --os ${output}/${base1}_${base2}.gen_samples

#calculate relatedness from sequence file
# vcftools --gzvcf ${file} --relatedness2 --out ${file}.pop_info

#gzip files
# gzip ${file}

#compare chromosomes
# bcftools stats /lustre/scratch113/projects/esgi-vbseq/20140319/20140320_VQSR2.5_reapply_9764/20140322_SHAPEIT/${file}.vcf.gz /lustre/scratch113/projects/fvg_seq/20140319/20140321_VQSR2.5_reapply_9654/20140323_SHAPEIT/${file}.vcf.gz > INGI_VB_FVG_comp_chr${file}.log

# plot stats
# plot-vcfstats INGI_VB_FVG_comp_chr${file}.log -p PLOTS/plot_${file}

# - Comparison between INGI and UK10K
# - FVG
# bcftools stats /lustre/scratch113/projects/fvg_seq/20140319/20140321_VQSR2.5_reapply_9654/20140323_SHAPEIT/${file}.vcf.gz /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/vcf_sites_filtered/chr${file}.sites.vcf.gz > INGI_FVG_UK10K_comp_chr${file}.log
# plot-vcfstats INGI_FVG_UK10K_comp_chr${file}.log -p PLOTS/plot_${file}
	 
# - VBI
# bcftools stats /lustre/scratch113/projects/esgi-vbseq/20140319/20140320_VQSR2.5_reapply_9764/20140322_SHAPEIT/${file}.vcf.gz /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/vcf_sites_filtered/chr${file}.sites.vcf.gz > INGI_VB_UK10K_comp_chr${file}.log
# plot-vcfstats INGI_VB_UK10K_comp_chr${file}.log -p PLOTS/plot_${file}

# - Comparison between INGI and TGP
# - FVG
# bcftools stats /lustre/scratch113/projects/fvg_seq/20140319/20140321_VQSR2.5_reapply_9654/20140323_SHAPEIT/${file}.vcf.gz /lustre/scratch111/resources/variation/Homo_sapiens/grch37/1K_phase1_release_v3.20101123/ALL.chr${file}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz > /lustre/scratch113/teams/soranzo/users/mc14/INGI/SEQ_COMPARISON_TGP/FVG/INGI_FVG_TGP_comp_chr${file}.log
# plot-vcfstats /lustre/scratch113/teams/soranzo/users/mc14/INGI/SEQ_COMPARISON_TGP/FVG/INGI_FVG_TGP_comp_chr${file}.log -p /lustre/scratch113/teams/soranzo/users/mc14/INGI/SEQ_COMPARISON_TGP/FVG/PLOTS/plot_${file}
 
# # - VBI
# bcftools stats /lustre/scratch113/projects/esgi-vbseq/20140319/20140320_VQSR2.5_reapply_9764/20140322_SHAPEIT/${file}.vcf.gz /lustre/scratch111/resources/variation/Homo_sapiens/grch37/1K_phase1_release_v3.20101123/ALL.chr${file}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz > /lustre/scratch113/teams/soranzo/users/mc14/INGI/SEQ_COMPARISON_TGP/VBI/INGI_VB_TGP_comp_chr${file}.log
# plot-vcfstats /lustre/scratch113/teams/soranzo/users/mc14/INGI/SEQ_COMPARISON_TGP/VBI/INGI_VB_TGP_comp_chr${file}.log -p /lustre/scratch113/teams/soranzo/users/mc14/INGI/SEQ_COMPARISON_TGP/VBI/PLOTS/plot_${file}

# # - INGI vs ALL
# vcf-compare /lustre/scratch113/projects/esgi-vbseq/20140319/20140320_VQSR2.5_reapply_9764/20140322_SHAPEIT/${file}.vcf.gz /lustre/scratch113/projects/fvg_seq/20140319/20140321_VQSR2.5_reapply_9654/20140323_SHAPEIT/${file}.vcf.gz /lustre/scratch111/resources/variation/Homo_sapiens/grch37/1K_phase1_release_v3.20101123/ALL.chr${file}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/vcf_sites_filtered/chr${file}.sites.vcf.gz > /lustre/scratch113/teams/soranzo/users/mc14/INGI/SEQ_COMPARISON_ALL/INGI_vs_all_comp_chr${file}.log
# plot-vcfstats /lustre/scratch113/teams/soranzo/users/mc14/INGI/SEQ_COMPARISON_ALL/INGI_vs_all_comp_chr${file}.log -p /lustre/scratch113/teams/soranzo/users/mc14/INGI/SEQ_COMPARISON_ALL/PLOTS/plot_${file}

#extract unrelated samples in the same way as we did for village split
# $2=vcf file path
# $3=outpath
# village_file=`basename ${file}`
# village=${village_file%_*}


# bcftools view -s ${file} $2 -O v | vcf-annotate --fill-ICF | bgzip -c > $2.${village}.vcf.gz
# tabix -f -p vcf $2.${village}.vcf.gz

# # extract stats by village
# bcftools stats -d 0,5000,1 -s - $2.${village}.vcf.gz > $2.${village}.vcfchk;plot-vcfstats $2.${village}.vcfchk -p vcf_check/unrelated/plots/${village}

# # extract table info by village
# #removed multiallelic sites and formatted the output to avoid missing problem when checking frequencies
# (echo "CHROM POS ID REF ALT AC AN TGP_AF TGP_AMR_AF TGP_ASN_AF TGP_AFR_AF TGP_EUR_AF IMP2 VQSLOD AF MAF MINOR";bcftools2 query $2.${village}.vcf.gz -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/1kg_AF\t%INFO/1kg_AMR_AF\t%INFO/1kg_ASN_AF\t%INFO/1kg_AFR_AF\t%INFO/1kg_EUR_AF\t%INFO/IMP2\t%INFO/VQSLOD\n'| fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " " > $2.${village}.vcf.gz.maf_table.tab

#extract stats after annotation
# bcftools stats -i ${file}.vcf.gz > ${file}.novel.vchk

#extract info using samtools on coverage
# filename=`basename ${file}`
# samtools mpileup -D  ${file} | awk '{x+=$4;next}END{print x/NR}' > ${filename}.mean_coverage
#use a higher depth cut off value to include all variants
# samtools mpileup -D -BQ0 -d10000000 ${file} | awk '{x+=$4;next}END{print x/NR}' > ${filename}.mean_coverage


#extract stats from TGP files
# filename=`basename ${file}`
# bcftools stats ${file} > ${filename}.novel.vchk

#extract table info by village
#removed multiallelic sites and formatted the output to avoid missing problem when checking frequencies
# (echo "CHROM POS ID REF ALT AC AN IMP2 VQSLOD AF MAF MINOR";bcftools2 query ${file}.vcf.gz -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/IMP2\t%INFO/VQSLOD\n'| fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " " > MAF/${file}.maf_table.tab


#convert SHAPEIT phased haplotypes to vcf for VBI population to use in NRD calculations
# shapeit2 -convert --input-haps /lustre/scratch113/projects/uk10k/users/jh21/imputed/vb/shapeit/chr${file}.hap.gz /lustre/scratch113/projects/uk10k/users/jh21/imputed/vb/shapeit/vbi.sample --output-vcf /lustre/scratch113/teams/soranzo/users/mc14/INGI_VB/ARRAY/chr${file}.vcf
# (grep "^#" /lustre/scratch113/teams/soranzo/users/mc14/INGI_VB/ARRAY/chr${file}.vcf;grep -v "^#" /lustre/scratch113/teams/soranzo/users/mc14/INGI_VB/ARRAY/chr${file}.vcf | awk -v chr=${file} '{print chr,$0}'|tr "\t" " "|cut -f 1,3- -d " "|tr " " "\t" ) | vcf-annotate --fill-ICF | bgzip -c > /lustre/scratch113/teams/soranzo/users/mc14/INGI_VB/ARRAY/chr${file}.vcf.gz
# tabix -p vcf /lustre/scratch113/teams/soranzo/users/mc14/INGI_VB/ARRAY/chr${file}.vcf.gz

#annotate using VEP fromcommand line. No beautyfy pipeline!!!
# /nfs/team151/software/variant_effect_predictor/variant_effect_predictor.pl -i ${file} --quiet --regulatory --sift b --polyphen b --plugin Condel,/software/vertres/bin-external/VEP_plugins/config/Condel/config/,b --symbol --format vcf --force_overwrite --cache --dir /nfs/team151/software/VEP

# extract table info and write table in current folder
# filename=`basename ${file}`
# (echo "CHROM POS ID REF ALT AN AC TGP_AF AMR_AF ASN_AF AFR_AF EUR_AF IMP2 VQSLOD TGP UK10K HWE DP ICF AF MAF MINOR";bcftools query ${file} -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\t%INFO/1kg_AF\t%INFO/1kg_AMR_AF\t%INFO/1kg_ASN_AF\t%INFO/1kg_AFR_AF\t%INFO/1kg_EUR_AF\t%INFO/IMP2\t%INFO/VQSLOD\t%INFO/TGP\t%INFO/UK10K\t%INFO/HWE\t%INFO/DP\t%INFO/ICF\n" | awk '{if($5 !~ ",") print $0}' | awk '{if($6 != 0)print $0,$7/$6;else print $0,"NA"}' | awk '{if($NF != "NA") {if($NF > 0.5) print $0,1-$NF,$4;else print $0,$NF,$5}else{print $0,"NA","NA"}}') | tr "\t" " " > ${filename}.csv
# (echo "CHROM POS ID REF ALT AN AC VQSLOD CULPRIT AF MAF MINOR";bcftools2 query ${file} -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\t%INFO/VQSLOD\t%INFO/culprit\n" | awk '{if($5 !~ ",") print $0}' | awk '{if($6 != 0)print $0,$7/$6;else print $0,"NA"}' | awk '{if($NF != "NA") {if($NF > 0.5) print $0,1-$NF,$4;else print $0,$NF,$5}else{print $0,"NA","NA"}}') | tr "\t" " " > ${file}.csv

#22/05/2014
#extract stats for vcf files
# filename=`basename ${file}`
# bcftools stats ${file} > ${filename}.vchk


# (echo "CHROM POS ID REF ALT AN AC TGP_AF AMR_AF ASN_AF AFR_AF EUR_AF IMP2 VQSLOD TGP UK10K HWE DP ICF AF MAF MINOR";bcftools query X.vcf.gz -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\t%INFO/1kg_AF\t%INFO/1kg_AMR_AF\t%INFO/1kg_ASN_AF\t%INFO/1kg_AFR_AF\t%INFO/1kg_EUR_AF\t%INFO/IMP2\t%INFO/VQSLOD\t%INFO/TGP\t%INFO/UK10K\t%INFO/HWE\t%INFO/DP\t%INFO/ICF\n" | awk '{if($5 !~ ",") print $0}' | awk '{if($6 != 0)print $0,$7/$6;else print $0,"NA"}' | awk '{if($NF != "NA") {if($NF > 0.5) print $0,1-$NF,$4;else print $0,$NF,$5}else{print $0,"NA","NA"}}') | tr "\t" " " > X.vcf.gz.csv

# extract table info and write table in current folder
# filename=`basename ${file}`
# (echo "CHROM POS ID REF ALT AN AC AF MAF MINOR";bcftools query ${file} -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\n" | awk '{if($5 !~ ",") print $0}' | awk '{if($6 != 0)print $0,$7/$6;else print $0,"NA"}' | awk '{if($NF != "NA") {if($NF > 0.5) print $0,1-$NF,$4;else print $0,$NF,$5}else{print $0,"NA","NA"}}') | tr "\t" " " > ${filename}.csv
# (echo "CHROM POS ID REF ALT AN AC VBI FVG AF MAF MINOR";bcftools query ${file} -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\t%INFO/VBI\t%INFO/FVG\n" | awk '{if($5 !~ ",") print $0}' | awk '{if($6 != 0)print $0,$7/$6;else print $0,"NA"}' | awk '{if($NF != "NA") {if($NF > 0.5) print $0,1-$NF,$4;else print $0,$NF,$5}else{print $0,"NA","NA"}}') | tr "\t" " " > ${filename}.csv
# (echo "CHROM POS ID REF ALT AN AC VQSLOD CULPRIT AF MAF MINOR";bcftools2 query ${file} -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\t%INFO/VQSLOD\t%INFO/culprit\n" | awk '{if($5 !~ ",") print $0}' | awk '{if($6 != 0)print $0,$7/$6;else print $0,"NA"}' | awk '{if($NF != "NA") {if($NF > 0.5) print $0,1-$NF,$4;else print $0,$NF,$5}else{print $0,"NA","NA"}}') | tr "\t" " " > ${file}.csv

# #5/06/2014
# filename=`basename ${file}`
# #Extract plink file format from wgs from TSI and CEU using a keplist
# plink2 --vcf ${file} --keep $2 --make-bed --keep-allele-order --double-id --biallelic-only list --snps-only --out ${filename}
# awk '{if($2==".") print $1,"chr"$1":"$4,$3,$4,$5,$6;else print $0}' ${filename}.bim | tr " " "\t" > ${filename}.bim.sanitized
# mv ${filename}.bim ${filename}.bim.old
# mv ${filename}.bim.sanitized ${filename}.bim


#2/07/2014 extract data from 1000G phase1_release_v3 for different populations splitted by chromosomes
# pop=$2

# bcftools view -S /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/${pop}/${pop}.keeplist -O z -o /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/${pop}/${file}.vcf.gz /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/WGS/ALL.chr${file}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz

#25/07/2014
# Join chromosomes from different populations

# mv /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/CEU/${file}.vcf.gz /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/CEU/${file}.vcf.old.gz;zcat /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/CEU/${file}.vcf.old.gz | sed "s/ID=GL,Number=./ID=GL,Number=G/"|sed "s/ID=AC,Number=./ID=AC,Number=A/" | bgzip -c > /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/CEU/${file}.vcf.gz;
# mv /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/TSI/${file}.vcf.gz /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/TSI/${file}.vcf.old.gz;zcat /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/TSI/${file}.vcf.old.gz | sed "s/ID=GL,Number=./ID=GL,Number=G/"|sed "s/ID=AC,Number=./ID=AC,Number=A/" | bgzip -c > /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/TSI/${file}.vcf.gz;

# tabix -f -p vcf /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/CEU/${file}.vcf.gz
# tabix -f -p vcf /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/TSI/${file}.vcf.gz

# bcftools merge /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/CEU/${file}.vcf.gz /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/TSI/${file}.vcf.gz /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/INGI/FVG/${file}.vcf.gz /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/INGI/VBI/${file}.vcf.gz /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/INGI/CARL/${file}.vcf.gz -O z > chr${file}.union.vcf.gz

#28/07/2014
#index joint files
# tabix -p vcf chr${file}.union.vcf.gz

# bcftools merge /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/CEU/21.vcf.gz /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/TGP/TSI/21.vcf.gz /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/INGI/FVG/21.vcf.gz /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/INGI/VBI/21.vcf.gz /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POPULATIONS/INGI/CARL/21.vcf.gz


#06/04/2014
#create bed files from vcf
# chr=${file}

#create output folder
# mkdir -p CHR${chr}

# ec_dacmacdafmaf2bed.py /lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/${pop}_unrelated.list /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140711_ANNOTATED/${chr}.vcf.gz ${pop}.chr${chr}.tab

# mv ${pop}.chr${chr}.tab CHR${chr}/
# gzip CHR${chr}/${pop}.chr${chr}.tab

#19/08/2014
#Create plink input files for NON MISSING data
# chr=${file}
# pop=$2

# if [ ! -s /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/ALL.${chr}.non_missing.ped ]
# then

#   plink2 --vcf /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING/${chr}.non_missing.vcf.gz --cm-map /nfs/team151/reference/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr${chr}_combined_b37.txt ${chr} --maf 0.01 --double-id --biallelic-only --snps-only --keep-allele-order --recode --out /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/${pop}.${chr}.non_missing

# fi
#22/08/2014
# remove genotypes from all population to use funseq for annotation
# bcftools view -G -O v -o ${file}.clean_annotated.nogeno.vcf /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POP_MERGED_FILES/FIVE_POPS/20140730_ANNOTATED/${file}.clean_annotated.vcf.gz

#annotate using funseq
# funseq.sh -f /lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POP_MERGED_FILES/FIVE_POPS/20140822_NO_GENO/${file}.clean_annotated.nogeno.vcf -m 1 -maf 0 -inf vcf -outf bed -nc

#create a table for annotation of the complete files
# (head -1 out/${file}.FunSEQ.bed|awk '{print toupper($0)}'|tr ";" "\t";tail -n+2 out/${file}.FunSEQ.bed|tr ";" "\t") | bgzip -c > out/${file}.FunSEQ.bed.gz

# 9/09/2014
# remove genotypes from all population to use funseq for annotation
# chr=${file}
# p_file=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/${pop}_unrelated.list

# bcftools view -G -O z -o ${pop}.${chr}.nogeno.vcf.gz -S ${p_file} /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140711_ANNOTATED/${chr}.vcf.gz

# 12/09/2014
# chr=${file}
# pop=$2

# if [ ! -s /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/FILTERED/ALL.${chr}.non_missing.filtered.ped ]
# then

#   plink2 --file /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/ALL.${chr}.non_missing --extract /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/IBD/GERMLINE/ALL_TOGETHER/CHR${chr}/ALL.${chr}.non_missing.match.keepsnps --cm-map /nfs/team151/reference/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr${chr}_combined_b37.txt ${chr} --maf 0.01 --double-id --biallelic-only --snps-only --keep-allele-order --recode --out /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/FILTERED/${pop}.${chr}.non_missing.filtered

# fi

# 12/09/2014
# create files of filtered regions in splitted files
# g++ hom_to_ped.cpp -o hom_to_ped
# chr=${file}
# MATCH=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/IBD/GERMLINE/ALL_TOGETHER/CHR${chr}/ALL.${chr}.non_missing.match
# # pop=$2
# for reg_file in `ls ${MATCH}.shareDens_R*.to_include.keepsnps`
# do
#   rs_start=`cut -f 1 -d " " ${reg_file}`
#   rs_end=`cut -f 2 -d " " ${reg_file}`

#   cat /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/ALL.${chr}.non_missing.ped | ped_to_hom > /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/PED2HOM/ALL.${chr}.non_missing.hom.ped
#   cp /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/ALL.${chr}.non_missing.map /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/PED2HOM/ALL.${chr}.non_missing.hom.map

#   plink2 --file /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/PED2HOM/ALL.${chr}.non_missing.hom --from ${rs_start} --to ${rs_end} --recode 12 --out /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/PED2HOM/ALL.${chr}.non_missing_${rs_start}-${rs_end}.hom

#   # g++ ped_to_hom.cpp -o ped_to_hom
#   cat /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/PED2HOM/ALL.${chr}.non_missing_${rs_start}-${rs_end}.hom.ped | hom_to_ped > /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/HOM2PED/ALL.${chr}.non_missing_${rs_start}-${rs_end}.ped
#   mv /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/PED2HOM/ALL.${chr}.non_missing_${rs_start}-${rs_end}.hom.map /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/HOM2PED/ALL.${chr}.non_missing_${rs_start}-${rs_end}.map
# done

#17/09/2014
#update cm map position
# chr=${file}

# if [ ! -s /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/MAP_UPDATED/ALL.${chr}.non_missing.ped ]
# then

# plink2 --file /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/ALL.${chr}.non_missing --cm-map /nfs/team151/reference/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr${chr}_combined_b37.txt ${chr} --maf 0.01 --double-id --biallelic-only --snps-only --keep-allele-order --recode --out /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/MAP_UPDATED/ALL.${chr}.non_missing

# fi

#18/09/2014
#assign pairs of IBD match to different population files
region=${file}
pop_file=$2
# pop_file=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/ALL_unrelated_8cohort.list
# pop_file=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/ALL_unrelated_5cohort.list
awk -v cohort=${pop_file} '
BEGIN{
  while (getline < cohort) {
    pop[$1]=$2;
  }
}
{ print $0,pop[$1]"_"pop[$3] }' ${region} | tr " " "\t" > ${region}.pop_added
# { print $0,pop[$1]"_"pop[$3] }' ALL.9.non_missing_rs75838906-rs9792573.ped.filtered.match | tr " " "\t" > ALL.9.non_missing_rs75838906-rs9792573.ped.filtered.match.pop_added
# for group in `awk '{print $(NF)}' ALL.9.non_missing_rs75838906-rs9792573.ped.filtered.match.pop_added|sort| uniq`
region_name=`basename ${region}`
for group in `awk '{print $(NF)}' ${region}.pop_added |sort| uniq`
do 
  fgrep -w ${group} ${region}.pop_added > /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/IBD/GERMLINE/FILTERED/POP_SPLIT/${region_name}.${group}
done

#clean up a littlebit
rm ${region}.pop_added
