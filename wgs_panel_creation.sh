#!/usr/local/bin/bash
#
# Script to select variants for reference panel creation
# it works with multiple populations and by chromosomes' vcf files, to parallelize all
# test

# vcf=/lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/22.vcf.gz
# cohort="VBI"
# outdir=/lustre/scratch113/projects/esgi-vbseq/27112015_INGI_REF_PANEL/
#ARGS:
# $1= vcf file for a single chromosome, better if named as "[chr].vcf.gz"
# $2=cohort
# $3=output folder
# 
# /nfs/users/nfs_m/mc14/Work/bash_scripts/wgs_panel_creation.sh /lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/22.vcf.gz VBI /lustre/scratch113/projects/esgi-vbseq/27112015_INGI_REF_PANEL/
vcf=$1
outdir=$3/$2

mkdir -p ${outdir}
filename=`basename ${vcf}`
first_suffix="${filename%%.*}"
#we're going to split snps and indels, than put them back together again

########## SNPS
#select sites with AC >=2 and DP>=5
bsub -J"extract_snps_acgt2dpgt5_${first_suffix}" -o"%J_extract_snps_acgt2dpgt5_${first_suffix}.o" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal -- bcftools query -i'TYPE="snp" && AC>=2 && DP>=5' -f "%CHROM\t%POS\n" ${vcf} -o ${outdir}/${filename}.snps_ac2dp5.tab
bsub -J"extract_snps_aceq0dpgt5_${first_suffix}" -o"%J_extract_snps_aceq0dpgt5_${first_suffix}.o" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal -- bcftools query -i'TYPE="snp" && AC==0 && DP>=5' -f "%CHROM\t%POS\n" ${vcf} -o ${outdir}/${filename}.snps_ac0dp5.tab

#select sites with AC=1 and DP>5
bsub -J"extract_snps_aceq1dpgt5_${first_suffix}" -o"%J_extract_snps_aceq1dpgt5_${first_suffix}.o" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal -- bcftools query -i"TYPE='snp' && AC==1 && DP>=5" -f "%CHROM\t%POS\n" ${vcf} -o ${outdir}/${filename}.snps_ac1dp5.tab

#test:
#VBI on chr 22 : removed 2 snps for DP<5

#extract a list of sites with AC1 in common between at least two isolates
#we already have the isec files for INGI here: /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/INGI/UNION/
#we need to check if the sum in col 5 of the sites.txt file, is >=2 (it means that my variant is shared at least between 2 isolates)
# cohort=sys.argv[1]
# var_list=sys.argv[2]
# overlap_list=sys.argv[3]
# outdir=sys.argv[4]
chr=${first_suffix}

bsub -J"extract_common_ingi_snps_aceq1dpgt5_${first_suffix}_INGI" -o"%J_extract_common_ingi_snps_aceq1dpgt5_${first_suffix}_INGI.o" -w "ended(extract_snps_aceq1dpgt5_${first_suffix})" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal --  python /nfs/users/nfs_m/mc14/Work/bash_scripts/overlap_check.py INGI ${outdir}/${filename}.snps_ac1dp5.tab /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/INGI/UNION/${chr}/sites.txt ${outdir}
bsub -J"extract_common_ingi_snps_aceq1dpgt5_${first_suffix}_TGP3" -o"%J_extract_common_ingi_snps_aceq1dpgt5_${first_suffix}_TGP3.o" -w "ended(extract_common_ingi_snps_aceq1dpgt5_${first_suffix}_INGI)" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal --  python /nfs/users/nfs_m/mc14/Work/bash_scripts/overlap_check.py TGP3 ${outdir}/${filename}.snps_ac1dp5.tab.not_over.tab /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/TGPph3/UNION/${chr}/sites.txt ${outdir}
bsub -J"extract_common_ingi_snps_aceq1dpgt5_${first_suffix}_UK10K" -o"%J_extract_common_ingi_snps_aceq1dpgt5_${first_suffix}_UK10K.o" -w "ended(extract_common_ingi_snps_aceq1dpgt5_${first_suffix}_INGI)" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal --  python /nfs/users/nfs_m/mc14/Work/bash_scripts/overlap_check.py UK10K ${outdir}/${filename}.snps_ac1dp5.tab.not_over.tab /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/UK10K/UNION/${chr}/sites.txt ${outdir}


#extract a list of sites in common with 1000G or/and UK10K


#merge those lists in a list of uniqe sites


#extract sites to retain



########## INDELS 
#select sites with AC >=2 and DP>=5
bsub -J"extract_indels" -o"%J_extract_indels_acgt2dpgt5.o" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal -- bcftools query -i'TYPE="indels" && AC>=2 && DP>=5' -f "%CHROM\t%POS\n" ${vcf} -o ${outdir}/${filename}.indels_ac2dp5.tab
bsub -J"extract_indels" -o"%J_extract_indels_aceq0dpgt5.o" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal -- bcftools query -i'TYPE="indels" && AC==0 && DP>=5' -f "%CHROM\t%POS\n" ${vcf} -o ${outdir}/${filename}.indels_ac0dp5.tab

#select sites with AC=1 and DP>5
bsub -J"extract_indels" -o"%J_extract_indels_aceq1dpgt5.o" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" -q normal -- bcftools query -i"TYPE='indels' && AC==1 && DP>=5" -f "%CHROM\t%POS\n" ${vcf} -o ${outdir}/${filename}.indels_ac1dp5.tab

#extract a list of sites with AC1 in common between at least two isolates

#extract a list of sites in common with 1000G or/and UK10K

#merge those lists in a list of uniqe sites

#extract sites to retain


