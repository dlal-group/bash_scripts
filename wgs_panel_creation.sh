#!/usr/bin/env bash
#
# Script to select variants for reference panel creation
# it works with multiple populations and by chromosomes' vcf files, to parallelize all
# test

# vcf=/lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/22.vcf.gz
# cohort="VBI"
# outdir=/lustre/scratch113/projects/esgi-vbseq/27112015_INGI_REF_PANEL/VBI
#ARGS:
# $1= vcf file for a single chromosome, better if named as "[chr].vcf.gz"
# $2=cohort
# $3=output folder
# $4=mode (snps/indels)

# /nfs/users/nfs_m/mc14/Work/bash_scripts/wgs_panel_creation.sh /lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/22.vcf.gz VBI /lustre/scratch113/projects/esgi-vbseq/27112015_INGI_REF_PANEL
set -e

vcf=$1
cohort=$2
outdir=$3/$2
mode=$4
stage=$5

mkdir -p ${outdir}/PANEL
mkdir -p ${outdir}/LOG_LISTS
mkdir -p ${outdir}/LOG_OVERLAPS
mkdir -p ${outdir}/LOG_PANEL

filename=`basename ${vcf}`
first_suffix="${filename%%.*}"
chr=${first_suffix}

mkdir -p ${outdir}/${chr}
#we're going to split snps and indels, than put them back together again

case ${stage} in
SELECT )
# define path
ingi_union=/lustre/scratch113/projects/esgi-vbseq/01032016_PANEL_SESOURCES/INGI/UNION
tgp_ingi=/lustre/scratch113/projects/esgi-vbseq/01032016_PANEL_SESOURCES/TGPph3/UNION
uk10k_ingi=/lustre/scratch113/projects/esgi-vbseq/01032016_PANEL_SESOURCES/UK10K/UNION
#select sites with AC >=2 and DP>=5
bsub -J"extract_${mode}_${cohort}_acgt2dpgt5_${first_suffix}" -o"${outdir}/LOG_LISTS/1.%J_extract_${mode}_${cohort}_acgt2dpgt5_${first_suffix}.o" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- bcftools query -i"TYPE='${mode}' && AC>=2 && DP>=5" -f "%CHROM\t%POS\t%REF\t%ALT\n" ${vcf} -o ${outdir}/${chr}/${filename}.${mode}_ac2dp5.tab
bsub -J"extract_${mode}_${cohort}_aceq0dpgt5_${first_suffix}" -o"${outdir}/LOG_LISTS/2.%J_extract_${mode}_${cohort}_aceq0dpgt5_${first_suffix}.o" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- bcftools query -i"TYPE='${mode}' && AC==0 && DP>=5" -f "%CHROM\t%POS\t%REF\t%ALT\n" ${vcf} -o ${outdir}/${chr}/${filename}.${mode}_ac0dp5.tab

#select sites with AC=1 and DP>5
bsub -J"extract_${mode}_${cohort}_aceq1dpgt5_${first_suffix}" -o"${outdir}/LOG_LISTS/3.%J_extract_${mode}_${cohort}_aceq1dpgt5_${first_suffix}.o" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- bcftools query -i"TYPE='${mode}' && AC==1 && DP>=5" -f "%CHROM\t%POS\t%REF\t%ALT\n" ${vcf} -o ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab

#test
#VBI on chr 22 : removed 2 snps for DP<5

#extract a list of sites with AC1 in common between at least two isolates
#we already have the isec files for INGI here: /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/INGI/UNION/
#we need to check if the sum in col 5 of the sites.txt file, is >=2 (it means that my variant is shared at least between 2 isolates)
# cohort=sys.argv[1]
# var_list=sys.argv[2]
# overlap_list=sys.argv[3]
# outdir=sys.argv[4]

bsub -J"extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_INGI" -o"${outdir}/LOG_OVERLAPS/4.%J_extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_INGI.o" -w "ended(extract_${mode}_${cohort}_aceq1dpgt5_${first_suffix})" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- python /nfs/users/nfs_m/mc14/Work/bash_scripts/overlap_check.py INGI ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab ${ingi_union}/${chr}/sites.txt ${outdir}/${chr} ${mode}

#extract a list of sites in common with 1000G or/and UK10K
bsub -J"extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_TGP3" -o"${outdir}/LOG_OVERLAPS/5.%J_extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_TGP3.o" -w "ended(extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_INGI)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- python /nfs/users/nfs_m/mc14/Work/bash_scripts/overlap_check.py TGP3 ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab ${tgp_ingi}/${chr}/sites.txt ${outdir}/${chr} ${mode}
bsub -J"extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_UK10K" -o"${outdir}/LOG_OVERLAPS/6.%J_extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_UK10K.o" -w "ended(extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_INGI)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- python /nfs/users/nfs_m/mc14/Work/bash_scripts/overlap_check.py UK10K ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab ${uk10k_ingi}/${chr}/sites.txt ${outdir}/${chr} ${mode}

#after all the pipeline run, we end up with the following files for each chr:
# 2864 9.vcf.gz.snp_ac0dp5.tab
# 28737 9.vcf.gz.snp_ac1dp5.tab.INGI.snp.over.tab
# 633910 9.vcf.gz.snp_ac2dp5.tab
# -> To keep as they are

# 157993 9.vcf.gz.snp_ac1dp5.tab
# to filter and split in the following ->

# 38099 9.vcf.gz.snp_ac1dp5.tab.INGI.snp.not_over.tab.TGP3.snp.over.tab
# 31237 9.vcf.gz.snp_ac1dp5.tab.INGI.snp.not_over.tab.UK10K.snp.over.tab
# To join together and keep uniqes records -> we get 48873 Uniq sites from TGP3 and UK10K

# 129256 9.vcf.gz.snp_ac1dp5.tab.INGI.snp.not_over.tab
# 91157 9.vcf.gz.snp_ac1dp5.tab.INGI.snp.not_over.tab.TGP3.snp.not_over.tab
# 98019 9.vcf.gz.snp_ac1dp5.tab.INGI.snp.not_over.tab.UK10K.snp.not_over.tab
# To throw away

# In  total we keep 714384 out of 857987 -> removing all singletons plus a few sites with DP < 5

#so we need to merge selected list files in a list of uniqe sites
# 9.vcf.gz.snp_ac0dp5.tab
# 9.vcf.gz.snp_ac2dp5.tab
# 9.vcf.gz.snp_ac1dp5.tab.INGI.snp.over.tab -> here we also have the information about the site being Biallelic or Multiallelic
# 9.vcf.gz.snp_ac1dp5.tab.INGI.snp.not_over.tab.TGP3.snp.over.tab -> here we also have the information about the site being Biallelic or Multiallelic
# 9.vcf.gz.snp_ac1dp5.tab.INGI.snp.not_over.tab.UK10K.snp.over.tab -> here we also have the information about the site being Biallelic or Multiallelic

#solution 1) extract sites to retain
#merge list for AC0 and AC2 and all fo the AC1 
# echo "cat ${outdir}/${chr}/${filename}.${mode}_ac0dp5.tab ${outdir}/${chr}/${filename}.${mode}_ac2dp5.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.over.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab.TGP3.${mode}.over.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab.UK10K.${mode}.over.tab | cut -f -4 | sort -g -k1,1 -k2,2 | uniq | cut -f -2 > ${outdir}/${chr}/${filename}.${mode}_extract_bed.list" | bsub -J"merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}" -o"${outdir}/LOG_OVERLAPS/7.%J_merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}.o" -w "ended(extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_*)" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal
# bsub -J"merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}" -o"%J_merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}.o" -w "ended(extract_common_ingi_${mode}_aceq1dpgt5_${first_suffix}_*)" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal -- (cat ${outdir}/${filename}.${mode}_ac0dp5.tab ${outdir}/${filename}.${mode}_ac2dp5.tab <(cat ${outdir}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.over.tab  ${outdir}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab.TGP3.${mode}.over.tab ${outdir}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab.UK10K.${mode}.over.tab | cut -f -4 )) | sort -g -k1,1 -k2,2 | uniq > ${outdir}/${filename}.${mode}_extract.list

# bsub -J"generate_ingi_${mode}_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/8.%J_generate_ingi_${mode}_${first_suffix}_${cohort}_vcf.o" -w"ended(merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort})" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal -- bcftools view -i"TYPE='${mode}'" -R ${outdir}/${chr}/${filename}.${mode}_extract.list -O z -o ${outdir}/PANEL/${filename}.${mode}_REF.vcf.gz ${vcf}
# bsub -J"generate_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi" -o"${outdir}/LOG_PANEL/8.%J_generate_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi.o" -w"ended(generate_ingi_${mode}_${first_suffix}_${cohort}_vcf)" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal -- tabix -p vcf ${outdir}/PANEL/${filename}.${mode}_REF.vcf.gz

;;
EXTRACT)
#solution 2) annotate sites to retain in original VCF than extract those sites only
echo "Sites extraction step!"
# echo "cat ${outdir}/${chr}/${filename}.${mode}_ac0dp5.tab ${outdir}/${chr}/${filename}.${mode}_ac2dp5.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.over.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab.TGP3.${mode}.over.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab.UK10K.${mode}.over.tab | cut -f -4 | sort -g -k1,1 -k2,2 | uniq| awk 'BEGIN{OFS=\"\t\"}{print \$0,\"REF_PANEL\"}' | bgzip -c > ${outdir}/${chr}/${filename}.${mode}_extract_annotation.list.gz" | bsub -J"merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}" -o"${outdir}/LOG_OVERLAPS/9.%J_merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}.o" -w "ended(extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_*)" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal
echo "cat ${outdir}/${chr}/${filename}.${mode}_ac0dp5.tab ${outdir}/${chr}/${filename}.${mode}_ac2dp5.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.over.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab.TGP3.${mode}.over.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab.UK10K.${mode}.over.tab | cut -f -4 | sort -g -k1,1 -k2,2 | uniq| awk 'BEGIN{OFS=\"\t\"}{print \$0,\"REF_PANEL\"}' | bgzip -c > ${outdir}/${chr}/${filename}.${mode}_extract_annotation.list.gz" | bsub -J"merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}" -o"${outdir}/LOG_OVERLAPS/9.%J_merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}.o" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal
bsub -J"merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}_tbi" -o"${outdir}/LOG_OVERLAPS/9.%J_merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}_tbi.o" -w "ended(merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort})" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- tabix -b 2 -e 2 -s 1 -f ${outdir}/${chr}/${filename}.${mode}_extract_annotation.list.gz

bsub -J"annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/10.%J_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf.o" -w"ended(merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}_tbi)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- bcftools annotate -a ${outdir}/${chr}/${filename}.${mode}_extract_annotation.list.gz -m REF_PANEL -c CHROM,POS,REF,ALT -O z -o ${outdir}/PANEL/${filename}.${mode}_annREF.vcf.gz ${vcf} 
bsub -J"annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi" -o"${outdir}/LOG_PANEL/10.%J_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi.o" -w"ended(annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- tabix -f -p vcf ${outdir}/PANEL/${filename}.${mode}_annREF.vcf.gz

bsub -J"extract_ingi_${mode}_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/11.%J_extract_ingi_${mode}_${first_suffix}_${cohort}_vcf.o" -w"ended(annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- bcftools view -i"TYPE='${mode}' && INFO/REF_PANEL=1" -O z -o ${outdir}/PANEL/${filename}.${mode}_REF.vcf.gz ${outdir}/PANEL/${filename}.${mode}_annREF.vcf.gz 
bsub -J"extract_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi" -o"${outdir}/LOG_PANEL/11.%J_extract_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi.o" -w"ended(extract_ingi_${mode}_${first_suffix}_${cohort}_vcf)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- tabix -f -p vcf ${outdir}/PANEL/${filename}.${mode}_REF.vcf.gz
;;
MERGE )
# add a bit for merg back SNP and INDEls
echo "Merge step!"
# bsub -J"merge_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/12.%J_merge_${first_suffix}_${cohort}_vcf.o" -w"ended(extract_ingi_snp_${first_suffix}_${cohort}_vcf_tbi) && ended(extract_ingi_indel_${first_suffix}_${cohort}_vcf_tbi)" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- bcftools concat -a  ${outdir}/PANEL/${filename}.snp_REF.vcf.gz ${outdir}/PANEL/${filename}.indel_REF.vcf.gz -O z -o ${outdir}/PANEL/${filename}.ALL_REF.vcf.gz
bsub -J"merge_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/12.%J_merge_${first_suffix}_${cohort}_vcf.o" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- bcftools concat -a  ${outdir}/PANEL/${filename}.snp_REF.vcf.gz ${outdir}/PANEL/${filename}.indel_REF.vcf.gz -O z -o ${outdir}/PANEL/${filename}.ALL_REF.vcf.gz
bsub -J"merge_${first_suffix}_${cohort}_vcf_tbi" -o"${outdir}/LOG_PANEL/12.%J_merge_${first_suffix}_${cohort}_vcf_tbi.o" -w"ended(merge_${first_suffix}_${cohort}_vcf)" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- tabix -p vcf ${outdir}/PANEL/${filename}.ALL_REF.vcf.gz
echo "bcftools stats -v -F /lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa ${outdir}/PANEL/${filename}.ALL_REF.vcf.gz > ${outdir}/PANEL/${filename}.ALL_REF.vcf.gz.vcfchk" | bsub -J"stats_merged_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/13.%J_stats_merged_${first_suffix}_${cohort}.o" -w"ended(merge_${first_suffix}_${cohort}_vcf_tbi)" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal

;;

CHECK )
# Retrieve duplicates line from extraction lists

# awk 'BEGIN{x_count[$1,$2]=1}
# {if (x[$1,$2]) {
#   x_count[$1,$2]++;
#    print $0;
#     if (x_count[$1,$2] == 2) {
#      print x[$1,$2];
#      print x_count[$1,$2]
#  }
#    }
#     x[$1,$2] = $0}'

# awk 'BEGIN{x_count[$1,$2]=0} {if (x[$1,$2]) {x_count[$1,$2]++; print x_count[$1,$2];if (x_count[$1,$2] == 2) { print x[$1,$2]; } } x[$1,$2] = $0}'


;;

ALL )

# define path
ingi_union=/lustre/scratch113/projects/esgi-vbseq/01032016_PANEL_SESOURCES/INGI/UNION
tgp_ingi=/lustre/scratch113/projects/esgi-vbseq/01032016_PANEL_SESOURCES/TGPph3/UNION
uk10k_ingi=/lustre/scratch113/projects/esgi-vbseq/01032016_PANEL_SESOURCES/UK10K/UNION
#select sites with AC >=2 and DP>=5
bsub -J"extract_${mode}_${cohort}_acgt2dpgt5_${first_suffix}" -o"${outdir}/LOG_LISTS/1.%J_extract_${mode}_${cohort}_acgt2dpgt5_${first_suffix}.o" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- bcftools query -i"TYPE='${mode}' && AC>=2 && DP>=5" -f "%CHROM\t%POS\t%REF\t%ALT\n" ${vcf} -o ${outdir}/${chr}/${filename}.${mode}_ac2dp5.tab
bsub -J"extract_${mode}_${cohort}_aceq0dpgt5_${first_suffix}" -o"${outdir}/LOG_LISTS/2.%J_extract_${mode}_${cohort}_aceq0dpgt5_${first_suffix}.o" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- bcftools query -i"TYPE='${mode}' && AC==0 && DP>=5" -f "%CHROM\t%POS\t%REF\t%ALT\n" ${vcf} -o ${outdir}/${chr}/${filename}.${mode}_ac0dp5.tab

#select sites with AC=1 and DP>5
bsub -J"extract_${mode}_${cohort}_aceq1dpgt5_${first_suffix}" -o"${outdir}/LOG_LISTS/3.%J_extract_${mode}_${cohort}_aceq1dpgt5_${first_suffix}.o" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- bcftools query -i"TYPE='${mode}' && AC==1 && DP>=5" -f "%CHROM\t%POS\t%REF\t%ALT\n" ${vcf} -o ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab

#test
#VBI on chr 22 : removed 2 snps for DP<5

#extract a list of sites with AC1 in common between at least two isolates
#we already have the isec files for INGI here: /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/INGI/UNION/
#we need to check if the sum in col 5 of the sites.txt file, is >=2 (it means that my variant is shared at least between 2 isolates)

bsub -J"extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_INGI" -o"${outdir}/LOG_OVERLAPS/4.%J_extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_INGI.o" -w "ended(extract_${mode}_${cohort}_aceq1dpgt5_${first_suffix})" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- python /nfs/users/nfs_m/mc14/Work/bash_scripts/overlap_check.py INGI ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab ${ingi_union}/${chr}/sites.txt ${outdir}/${chr} ${mode}

#extract a list of sites in common with 1000G or/and UK10K
bsub -J"extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_TGP3" -o"${outdir}/LOG_OVERLAPS/5.%J_extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_TGP3.o" -w "ended(extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_INGI)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- python /nfs/users/nfs_m/mc14/Work/bash_scripts/overlap_check.py TGP3 ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab ${tgp_ingi}/${chr}/sites.txt ${outdir}/${chr} ${mode}
bsub -J"extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_UK10K" -o"${outdir}/LOG_OVERLAPS/6.%J_extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_UK10K.o" -w "ended(extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_INGI)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- python /nfs/users/nfs_m/mc14/Work/bash_scripts/overlap_check.py UK10K ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab ${uk10k_ingi}/${chr}/sites.txt ${outdir}/${chr} ${mode}

#solution 2) annotate sites to retain in original VCF than extract those sites only
echo "Sites extraction step!"
echo "cat ${outdir}/${chr}/${filename}.${mode}_ac0dp5.tab ${outdir}/${chr}/${filename}.${mode}_ac2dp5.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.over.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab.TGP3.${mode}.over.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab.UK10K.${mode}.over.tab | cut -f -4 | sort -g -k1,1 -k2,2 | uniq| awk 'BEGIN{OFS=\"\t\"}{print \$0,\"REF_PANEL\"}' | bgzip -c > ${outdir}/${chr}/${filename}.${mode}_extract_annotation.list.gz" | bsub -J"merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}" -o"${outdir}/LOG_OVERLAPS/9.%J_merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}.o" -w "ended(extract_common_ingi_${mode}_${cohort}_aceq1dpgt5_${first_suffix}_*)" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal
# echo "cat ${outdir}/${chr}/${filename}.${mode}_ac0dp5.tab ${outdir}/${chr}/${filename}.${mode}_ac2dp5.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.over.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab.TGP3.${mode}.over.tab ${outdir}/${chr}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab.UK10K.${mode}.over.tab | cut -f -4 | sort -g -k1,1 -k2,2 | uniq| awk 'BEGIN{OFS=\"\t\"}{print \$0,\"REF_PANEL\"}' | bgzip -c > ${outdir}/${chr}/${filename}.${mode}_extract_annotation.list.gz" | bsub -J"merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}" -o"${outdir}/LOG_OVERLAPS/9.%J_merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}.o" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal
bsub -J"merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}_tbi" -o"${outdir}/LOG_OVERLAPS/9.%J_merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}_tbi.o" -w "ended(merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort})" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- tabix -b 2 -e 2 -s 1 -f ${outdir}/${chr}/${filename}.${mode}_extract_annotation.list.gz

bsub -J"annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/10.%J_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf.o" -w"ended(merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}_tbi)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- bcftools annotate -a ${outdir}/${chr}/${filename}.${mode}_extract_annotation.list.gz -m REF_PANEL -c CHROM,POS,REF,ALT -O z -o ${outdir}/PANEL/${filename}.${mode}_annREF.vcf.gz ${vcf} 
bsub -J"annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi" -o"${outdir}/LOG_PANEL/10.%J_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi.o" -w"ended(annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- tabix -f -p vcf ${outdir}/PANEL/${filename}.${mode}_annREF.vcf.gz

bsub -J"extract_ingi_${mode}_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/11.%J_extract_ingi_${mode}_${first_suffix}_${cohort}_vcf.o" -w"ended(annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- bcftools view -i"TYPE='${mode}' && INFO/REF_PANEL=1" -O z -o ${outdir}/PANEL/${filename}.${mode}_REF.vcf.gz ${outdir}/PANEL/${filename}.${mode}_annREF.vcf.gz 
bsub -J"extract_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi" -o"${outdir}/LOG_PANEL/11.%J_extract_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi.o" -w"ended(extract_ingi_${mode}_${first_suffix}_${cohort}_vcf)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- tabix -f -p vcf ${outdir}/PANEL/${filename}.${mode}_REF.vcf.gz

###########
echo "Cleaning step!!!"
#we need to remove all those normalized stuff that is confounding (that should result in removal of SOME indels sites, not snps)

#first create a keeping table in the form of
# CHROM POS REF ALT INFO/DP4 INFO/DP4 INFO/DP INFO/HOB INFO/ICB INFO/IDV 
bsub -J"clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist" -o"${outdir}/LOG_PANEL/13.%J_clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist.o" -w"ended(extract_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- bcftools query -f"%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP4\t%INFO/DP\t%INFO/HOB\t%INFO/ICB\t%INFO/IDV\n" ${outdir}/PANEL/${filename}.${mode}_REF.vcf.gz -o ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab
#now use the python script to mark sites to keep
bsub -J"clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_annotate" -o"${outdir}/LOG_PANEL/13.%J_clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_annotate.o" -w"ended(clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- python /nfs/users/nfs_m/mc14/Work/bash_scripts/panel_check.py ${cohort} ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab ${outdir}/PANEL/ ${mode}

#bgzip and index the thing
echo "sort -g -k2,2 ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab.${cohort}.${mode}.to_keep.tab| bgzip -c > ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab.${cohort}.${mode}.to_keep.tab.gz" | bsub -J"clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_gzip" -o"${outdir}/LOG_PANEL/13.%J_clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_gzip.o" -w"ended(clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_annotate)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal

bsub -J"clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_gzip_index" -o"${outdir}/LOG_PANEL/13.%J_clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_gzip_index.o" -w"ended(clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_gzip)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- tabix -f -s 1 -b 2 -e 2 ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab.${cohort}.${mode}.to_keep.tab.gz

#now use the indexed file to annotate the vcf and than extract the final variant set
bsub -J"clean_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/14.%J_clean_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf.o" -w"ended(clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_gzip_index)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- bcftools annotate -a ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab.${cohort}.${mode}.to_keep.tab.gz -m KEEP -c CHROM,POS,REF,ALT,INFO/DP4,INFO/DP,INFO/HOB,INFO/ICB,INFO/IDV -O z -o ${outdir}/PANEL/${filename}.${mode}_KEEP_REF.vcf.gz ${outdir}/PANEL/${filename}.${mode}_REF.vcf.gz
bsub -J"clean_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi" -o"${outdir}/LOG_PANEL/14.%J_clean_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi.o" -w"ended(clean_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- tabix -f -p vcf ${outdir}/PANEL/${filename}.${mode}_KEEP_REF.vcf.gz

bsub -J"extract_clean_ingi_${mode}_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/15.%J_extract_clean_ingi_${mode}_${first_suffix}_${cohort}_vcf.o" -w"ended(clean_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- bcftools view -i"TYPE='${mode}' && INFO/KEEP=1" -O z -o ${outdir}/PANEL/${filename}.${mode}_clean_REF.vcf.gz ${outdir}/PANEL/${filename}.${mode}_KEEP_REF.vcf.gz
bsub -J"extract_clean_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi" -o"${outdir}/LOG_PANEL/15.%J_extract_clean_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi.o" -w"ended(extract_clean_ingi_${mode}_${first_suffix}_${cohort}_vcf)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- tabix -f -p vcf ${outdir}/PANEL/${filename}.${mode}_clean_REF.vcf.gz


################################################
# add a bit for merg back SNP and INDEls
echo "Merge step!"
bsub -J"merge_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/12.%J_merge_${first_suffix}_${cohort}_vcf.o" -w"ended(extract_ingi_snp_${first_suffix}_${cohort}_vcf_tbi) && ended(extract_ingi_indel_${first_suffix}_${cohort}_vcf_tbi)" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- bcftools concat -a  ${outdir}/PANEL/${filename}.snp_REF.vcf.gz ${outdir}/PANEL/${filename}.indel_REF.vcf.gz -O z -o ${outdir}/PANEL/${filename}.ALL_REF.vcf.gz
# bsub -J"merge_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/12.%J_merge_${first_suffix}_${cohort}_vcf.o" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- bcftools concat -a  ${outdir}/PANEL/${filename}.snp_REF.vcf.gz ${outdir}/PANEL/${filename}.indel_REF.vcf.gz -O z -o ${outdir}/PANEL/${filename}.ALL_REF.vcf.gz
bsub -J"merge_${first_suffix}_${cohort}_vcf_tbi" -o"${outdir}/LOG_PANEL/12.%J_merge_${first_suffix}_${cohort}_vcf_tbi.o" -w"ended(merge_${first_suffix}_${cohort}_vcf)" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal -- tabix -p vcf ${outdir}/PANEL/${filename}.ALL_REF.vcf.gz
echo "bcftools stats -v -F /lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa ${outdir}/PANEL/${filename}.ALL_REF.vcf.gz > ${outdir}/PANEL/${filename}.ALL_REF.vcf.gz.vcfchk" | bsub -J"stats_merged_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/13.%J_stats_merged_${first_suffix}_${cohort}.o" -w"ended(merge_${first_suffix}_${cohort}_vcf_tbi)" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal

;;

CLEAN)
echo "Cleaning step!!!"
#we need to remove all those normalized stuff that is confounding (that should result in removal of SOME indels sites, not snps)

#first create a keeping table in the form of
# CHROM POS REF ALT INFO/DP4 INFO/DP4 INFO/DP INFO/HOB INFO/ICB INFO/IDV 
bsub -J"clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist" -o"${outdir}/LOG_PANEL/13.%J_clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist.o" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- bcftools query -f"%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP4\t%INFO/DP\t%INFO/HOB\t%INFO/ICB\t%INFO/IDV\n" ${outdir}/PANEL/${filename}.${mode}_REF.vcf.gz -o ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab
#now use the python script to mark sites to keep
bsub -J"clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_annotate" -o"${outdir}/LOG_PANEL/13.%J_clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_annotate.o" -w"ended(clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- python /nfs/users/nfs_m/mc14/Work/bash_scripts/panel_check.py ${cohort} ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab ${outdir}/PANEL/ ${mode}

#bgzip and index the thing
echo "sort -g -k2,2 ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab.${cohort}.${mode}.to_keep.tab| bgzip -c > ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab.${cohort}.${mode}.to_keep.tab.gz" | bsub -J"clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_gzip" -o"${outdir}/LOG_PANEL/13.%J_clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_gzip.o" -w"ended(clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_annotate)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal

bsub -J"clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_gzip_index" -o"${outdir}/LOG_PANEL/13.%J_clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_gzip_index.o" -w"ended(clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_gzip)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- tabix -f -s 1 -b 2 -e 2 ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab.${cohort}.${mode}.to_keep.tab.gz

#now use the indexed file to annotate the vcf and than extract the final variant set
bsub -J"clean_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/14.%J_clean_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf.o" -w"ended(clean_ingi_${mode}_${first_suffix}_${cohort}_keeplist_gzip_index)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- bcftools annotate -a ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab.${cohort}.${mode}.to_keep.tab.gz -m KEEP -c CHROM,POS,REF,ALT,INFO/DP4,INFO/DP,INFO/HOB,INFO/ICB,INFO/IDV -O z -o ${outdir}/PANEL/${filename}.${mode}_KEEP_REF.vcf.gz ${outdir}/PANEL/${filename}.${mode}_REF.vcf.gz
bsub -J"clean_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi" -o"${outdir}/LOG_PANEL/14.%J_clean_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi.o" -w"ended(clean_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- tabix -f -p vcf ${outdir}/PANEL/${filename}.${mode}_KEEP_REF.vcf.gz

bsub -J"extract_clean_ingi_${mode}_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG_PANEL/15.%J_extract_clean_ingi_${mode}_${first_suffix}_${cohort}_vcf.o" -w"ended(clean_annotate_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- bcftools view -i"TYPE='${mode}' && INFO/KEEP=1" -O z -o ${outdir}/PANEL/${filename}.${mode}_clean_REF.vcf.gz ${outdir}/PANEL/${filename}.${mode}_KEEP_REF.vcf.gz
bsub -J"extract_clean_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi" -o"${outdir}/LOG_PANEL/15.%J_extract_clean_ingi_${mode}_${first_suffix}_${cohort}_vcf_tbi.o" -w"ended(extract_clean_ingi_${mode}_${first_suffix}_${cohort}_vcf)" -M 2000 -R "select[mem>=2000] rusage[mem=2000]" -q normal -- tabix -f -p vcf ${outdir}/PANEL/${filename}.${mode}_clean_REF.vcf.gz

;;
esac
