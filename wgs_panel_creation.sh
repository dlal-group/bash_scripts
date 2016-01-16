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
vcf=$1
cohort=$2
outdir=$3/$2
mode=$4

mkdir -p ${outdir}/PANEL
mkdir -p ${outdir}/LOG
filename=`basename ${vcf}`
first_suffix="${filename%%.*}"
#we're going to split snps and indels, than put them back together again

########## SNPS
#select sites with AC >=2 and DP>=5
bsub -J"extract_${mode}_acgt2dpgt5_${first_suffix}" -o"${outdir}/LOG/%J_extract_${mode}_acgt2dpgt5_${first_suffix}.o" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal -- bcftools query -i"TYPE='${mode}' && AC>=2 && DP>=5" -f "%CHROM\t%POS\t%REF\t%ALT\n" ${vcf} -o ${outdir}/${filename}.${mode}_ac2dp5.tab
bsub -J"extract_${mode}_aceq0dpgt5_${first_suffix}" -o"${outdir}/LOG/%J_extract_${mode}_aceq0dpgt5_${first_suffix}.o" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal -- bcftools query -i"TYPE='${mode}' && AC==0 && DP>=5" -f "%CHROM\t%POS\t%REF\t%ALT\n" ${vcf} -o ${outdir}/${filename}.${mode}_ac0dp5.tab

#select sites with AC=1 and DP>5
bsub -J"extract_${mode}_aceq1dpgt5_${first_suffix}" -o"${outdir}/LOG/%J_extract_${mode}_aceq1dpgt5_${first_suffix}.o" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal -- bcftools query -i"TYPE='${mode}' && AC==1 && DP>=5" -f "%CHROM\t%POS\t%REF\t%ALT\n" ${vcf} -o ${outdir}/${filename}.${mode}_ac1dp5.tab

#test
#VBI on chr 22 : removed 2 snps for DP<5

#extract a list of sites with AC1 in common between at least two isolates
#we already have the isec files for INGI here: /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/INGI/UNION/
#we need to check if the sum in col 5 of the sites.txt file, is >=2 (it means that my variant is shared at least between 2 isolates)
# cohort=sys.argv[1]
# var_list=sys.argv[2]
# overlap_list=sys.argv[3]
# outdir=sys.argv[4]
chr=${first_suffix}

bsub -J"extract_common_ingi_${mode}_aceq1dpgt5_${first_suffix}_INGI" -o"${outdir}/LOG/%J_extract_common_ingi_${mode}_aceq1dpgt5_${first_suffix}_INGI.o" -w "ended(extract_${mode}_aceq1dpgt5_${first_suffix})" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal --  python /nfs/users/nfs_m/mc14/Work/bash_scripts/overlap_check.py INGI ${outdir}/${filename}.${mode}_ac1dp5.tab /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/INGI/UNION/${chr}/sites.txt ${outdir} ${mode}

#extract a list of sites in common with 1000G or/and UK10K
bsub -J"extract_common_ingi_${mode}_aceq1dpgt5_${first_suffix}_TGP3" -o"${outdir}/LOG/%J_extract_common_ingi_${mode}_aceq1dpgt5_${first_suffix}_TGP3.o" -w "ended(extract_common_ingi_${mode}_aceq1dpgt5_${first_suffix}_INGI)" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal --  python /nfs/users/nfs_m/mc14/Work/bash_scripts/overlap_check.py TGP3 ${outdir}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/TGPph3/UNION/${chr}/sites.txt ${outdir} ${mode}
bsub -J"extract_common_ingi_${mode}_aceq1dpgt5_${first_suffix}_UK10K" -o"${outdir}/LOG/%J_extract_common_ingi_${mode}_aceq1dpgt5_${first_suffix}_UK10K.o" -w "ended(extract_common_ingi_${mode}_aceq1dpgt5_${first_suffix}_INGI)" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal --  python /nfs/users/nfs_m/mc14/Work/bash_scripts/overlap_check.py UK10K ${outdir}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/UK10K/UNION/${chr}/sites.txt ${outdir} ${mode}

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

#merge list for AC0 and AC2 and all fo the AC1 
echo "cat \${outdir}/\${filename}.\${mode}_ac0dp5.tab \${outdir}/\${filename}.\${mode}_ac2dp5.tab \${outdir}/\${filename}.\${mode}_ac1dp5.tab.INGI.\${mode}.over.tab  \${outdir}/\${filename}.\${mode}_ac1dp5.tab.INGI.\${mode}.not_over.tab.TGP3.\${mode}.over.tab \${outdir}/\${filename}.\${mode}_ac1dp5.tab.INGI.\${mode}.not_over.tab.UK10K.\${mode}.over.tab | cut -f -4 | sort -g -k1,1 -k2,2 | uniq | cut -f -2 > \${outdir}/\${filename}.\${mode}_extract.list " | bsub -J"merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}" -o"${outdir}/LOG/%J_merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}.o" -w "ended(extract_common_ingi_${mode}_aceq1dpgt5_${first_suffix}_*)" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal
# bsub -J"merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}" -o"%J_merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort}.o" -w "ended(extract_common_ingi_${mode}_aceq1dpgt5_${first_suffix}_*)" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal -- (cat ${outdir}/${filename}.${mode}_ac0dp5.tab ${outdir}/${filename}.${mode}_ac2dp5.tab <(cat ${outdir}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.over.tab  ${outdir}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab.TGP3.${mode}.over.tab ${outdir}/${filename}.${mode}_ac1dp5.tab.INGI.${mode}.not_over.tab.UK10K.${mode}.over.tab | cut -f -4 )) | sort -g -k1,1 -k2,2 | uniq > ${outdir}/${filename}.${mode}_extract.list

#extract sites to retain
bsub -J"generate_ingi_${mode}_${first_suffix}_${cohort}_vcf" -o"${outdir}/LOG/%J_generate_ingi_${mode}_${first_suffix}_${cohort}_vcf.o" -w"ended(merge_common_ingi_${mode}_aceq0eq1gt2dpgt5_${first_suffix}_${cohort})" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal -- bcftools view -i"TYPE='${mode}'" -R ${outdir}/${filename}.${mode}_extract.list -O z -o ${outdir}/PANEL/${filename}.${mode}_REF.vcf.gz ${vcf}
tabix -p vcf ${outdir}/PANEL/${filename}.${mode}_REF.vcf.gz

#need to put back together SNP and INDELs for each chromosome after everything is done
