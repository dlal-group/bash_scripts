#!/usr/local/bin/bash
#script to check if data in the array are updated to the last position and on the correct strand for further comparisons
#ARGS:
bim_file_pref=$1
h3_bim_file_in=$2 #/nfs/team151/jh21/data/hapmap/hapmap.b37.bim
out_pref=$3
out2_pref=$4

plink --noweb --id-match ${bim_file_pref}.bim SNP,2 ${h3_bim_file_in} SNP,2 + complete --out ${out_pref}
awk '$1 != $7 {print $2, $7}' ${out_pref}.matched > ${out2_pref}.bad_chr ## NOT expected
awk '$4 != $10 {print $2, $10}' ${out_pref}.matched > ${out2_pref}.bad_pos
awk '($5$6 != $11$12) && ($5$6 != $12$11) {print $2}' ${out_pref}.matched > ${out2_pref}.bad_allele ## FLIP or EXCLUDE
# plink --noweb --bfile /nfs/team151/jh21/data/fvg_370/fvg_370 --update-map fvg_370.bad_chr --update-chr --make-bed --out fvg_370.good_chr
plink --noweb --bfile ${bim_file_pref} --update-map ${out2_pref}.bad_pos --flip ${out2_pref}.bad_allele --make-bed --out ${out2_pref}.flipped

plink --noweb --id-match ${out2_pref}.flipped.bim SNP,2 ${h3_bim_file_in} SNP,2 + complete --out ${out2_pref}.flipped_h3

awk '($5$6 != $11$12) && ($5$6 != $12$11) {print $2}' ${out2_pref}.flipped_h3.matched > ${out2_pref}.flipped.bad_allele
plink --noweb --bfile ${out2_pref}.flipped --exclude ${out2_pref}.flipped.bad_allele --make-bed --out ${out2_pref}

#we don't need all those check, not now, at least
# plink --noweb --bfile fvg_370 --check-sex --out fvg_370
# plink --noweb --bfile fvg_370 --missing --out fvg_370
# plink --noweb --bfile fvg_370 --exclude /nfs/team151/jh21/files/high-LD-regions.txt --range --indep-pairwise 50 05 0.25 --out fvg_370
# plink --noweb --bfile fvg_370 --extract fvg_370.prune.in --make-bed --out fvg_370.pruned
# plink --noweb --bfile fvg_370.pruned --genome --maf 0.1 --out fvg_370
# plink --noweb --bfile fvg_370.pruned --het --out fvg_370
# plink --noweb --bfile fvg_370.pruned --bmerge /nfs/team151/jh21/data/hapmap/hapmap3.b37.bed /nfs/team151/jh21/data/hapmap/hapmap3.b37.bim /nfs/team151/jh21/data/hapmap/hapmap3.b37.fam --geno 0.05 --make-bed --out ${out_pref}
# cp ${out_pref}.bim ${out_pref}.pedsnp
# awk '{print $1,$2,$3,$4,$5,1}' ${out_pref}.fam > ${out_pref}.pedind
# smartpca.perl -i ${out_pref}.bed -a ${out_pref}.pedsnp -b ${out_pref}.pedind -o ${out_pref}.pca -p ${out_pref}.plot -e ${out_pref}.eval -l ${out_pref}.pca.log -m 0 -k 4
# R CMD BATCH --vanilla --slave fvg_370.genotypeQC.R.CMD


#06/02/2017
# Snippets to check exome chip data on FVG and CARL after imputation
# 1) Get exome chip sites used in merged dataset for each cohort
for pop in FVG
do
for chr in {1..22}
do
awk '{OFS="\t"}{print $1,$4,$4}' ${chr}/chr${chr}_sorted.bim > /home/cocca/imputation/31012017_MERGED_TEST/${pop}/${chr}_exome_sites_table.txt
done
done

# 2) get frequencies from 1000G data, splitting multiallelic sites
for pop in FVG
do
for chr in {1..22}
do

bcftools view -v snps -R /home/cocca/imputation/31012017_MERGED_TEST/${pop}/${chr}_exome_sites_table.txt /netapp/nfs/resources/1000GP_phase3/vcf/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz |bcftools norm -m -  | bcftools query -f"%CHROM\t%POS\t%REF\t%ALT\t%ID\t%AN\t%AC\n" | awk '{OFS="\t"}{print $0,$7/$6}'| awk '{OFS="\t"}{if($(NF)<=0.5) print $0,$(NF);else print $0,1-$(NF)}'> /home/cocca/imputation/31012017_MERGED_TEST/${pop}/${chr}_exome_TGP_sites.txt

done
done

#run the script to extract frequencies and the formated file to perform the test on
for pop in CARL FVG
do
for chr in {1..22}
do
exome_chip="/netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${chr}/chr${chr}_sorted.bim"
tgp_data="/home/cocca/imputation/31012017_MERGED_TEST/${pop}/${chr}_exome_TGP_sites.txt"
impute_data="/home/cocca/imputation/MERGED_INGI_TGP3_23012017/${pop}/MERGED/ALL/chr${chr}.gen_info"
chromosome=${chr}
# pop="CARL"
outpath="/home/cocca/imputation/31012017_MERGED_TEST/${pop}"
echo "Processing $pop $chr..."
/home/cocca/scripts/bash_scripts/exome_sites_extraction.py ${exome_chip} ${tgp_data} ${impute_data} ${chr} ${pop} ${outpath}
done
done

#for each cohort remove the selected sites
for pop in CARL FVG
for pop in FVG
do
basefolder="/home/cocca/imputation/MERGED_INGI_TGP3_23012017"
mkdir -p ${basefolder}/${pop}/MERGED/ALL/UNFILTERED
mkdir -p ${basefolder}/${pop}/MERGED/ALL/FILTERED

for chr in {2..22}
do
echo "/home/cocca/scripts/bash_scripts/qctool_remove_script.sh ${pop} ${chr} ${basefolder}"|qsub -N ${pop}_${chr}_clean -o ${basefolder}/\$JOB_ID_${pop}_${chr}_clean.log -e ${basefolder}/\$JOB_ID_${pop}_${chr}_clean.e -V -l h_vmem=10G

done
done