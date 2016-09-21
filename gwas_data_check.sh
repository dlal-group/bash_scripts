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
