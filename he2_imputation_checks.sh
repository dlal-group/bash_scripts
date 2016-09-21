#!/usr/local/bin/bash

#PIPELINE SCIRPT FOR IMPUTATION CHECKs

# /lustre/scratch114/teams/soranzo/projects/Clean_Imputed_data/cbr_phase2
# /lustre/scratch114/teams/soranzo/projects/Clean_Imputed_data/cbr_phase2/plink_cbr2_chr1_recode.bed

# /lustre/scratch114/teams/soranzo/projects/Clean_Imputed_data/cbr_phase3

# /lustre/scratch113/teams/soranzo/users/he2/blood_GWAS_cleanup/cbr_phase4_imputation

# /lustre/scratch113/teams/soranzo/users/he2/blood_GWAS_cleanup/cbr_phase4_imputation/plink_cbr4_chr1_recode.bed

# #standard pruning with plink
# /nfs/team151/software/plink2_18_April_2015/plink --bfile /lustre/scratch114/teams/soranzo/projects/Clean_Imputed_data/cbr_phase2/plink_cbr2_chr1_recode --indep-pairwise 5000 1000 0.2 --out plink_cbr2_chr1_recode
# /nfs/team151/software/plink2_18_April_2015/plink --bfile /lustre/scratch114/teams/soranzo/projects/Clean_Imputed_data/cbr_phase2/plink_cbr2_chr1_recode --extract plink_cbr2_chr1_recode.prune.in --make-bed --out plink_cbr2_chr1_recode_pruned

# /nfs/team151/software/plink2_18_April_2015/plink --bfile /lustre/scratch113/teams/soranzo/users/he2/blood_GWAS_cleanup/cbr_phase4_imputation/plink_chr1_recode --indep-pairwise 5000 1000 0.2 --out plink_cbr4_chr1_recode
# /nfs/team151/software/plink2_18_April_2015/plink --bfile /lustre/scratch113/teams/soranzo/users/he2/blood_GWAS_cleanup/cbr_phase4_imputation/plink_chr1_recode --extract plink_cbr4_chr1_recode.prune.in --make-bed --out plink_cbr4_chr1_recode_pruned

#fix rs id problem on imputed data
# for CHR in {1..22}
# do

# echo "zcat chr\${CHR}.gen.gz | awk -v chr=\${CHR} '{print \"chr\"chr\":\"\$3,\$0}' | cut -f 1,4- -d \" \" | awk '{print \"---\",\$0}'  | gzip -c > chr\${CHR}.chrpos.gen.gz " | bsub -J"$CHR_chrpos" -o "%J_$CHR_chrpos.o" -q normal -R"select[mem>]
# echo "(head -1 chr\${CHR}.gen_info ;(fgrep -v position chr\${CHR}.gen_info | awk -v chr=\${CHR} '{print \"chr\"chr\":\"\$3,\$0}' | cut -f 1,4- -d \" \" | awk '{print \"---\",\$0}')) > chr\${CHR}.chrpos.gen_info " |  bsub -J"$CHR_chrpos" -o "%J_$CHR_chrpos.o" -q normal -R"select[mem>]

# done
# #fix rs id problem on imputed data plink format
# for CHR in {1..21}
# do
# echo $CHR
# awk '{print $1,"chr"$1":"$4,$3,$4,$5,$6}' cbr2_chr${CHR}_plink.bim | tr " " "\t" > cbr2_chr${CHR}_plink_mod.bim
# mv cbr2_chr${CHR}_plink.bim cbr2_chr${CHR}_plink_OLD.bim
# mv cbr2_chr${CHR}_plink_mod.bim cbr2_chr${CHR}_plink.bim
# done



CHR=1
outpath=/lustre/scratch113/teams/soranzo/users/mc14/HE2checks

# CHR=$1
# outpath=$2


#standard pruning on the fly with plink
cbr2_imp_path=/lustre/scratch114/teams/soranzo/projects/Clean_Imputed_data/cbr_phase2/imputed/uk10k1kg.shapeit
/nfs/team151/software/plink2_18_April_2015/plink --gen ${cbr2_imp_path}/chr${CHR}.gen.gz --sample ${cbr2_imp_path}/chr${CHR}.01.gen_samples --oxford-single-chr ${CHR} --indep-pairwise 5000 1000 0.1 --out ${outpath}/plink_cbr2_chr${CHR}_recode
/nfs/team151/software/plink2_18_April_2015/plink --gen ${cbr2_imp_path}/chr${CHR}.gen.gz --sample ${cbr2_imp_path}/chr${CHR}.01.gen_samples --oxford-single-chr ${CHR} --extract plink_cbr2_chr${CHR}_recode.prune.in --make-bed --out ${outpath}/plink_cbr2_chr${CHR}_recode_pruned
/nfs/team151/software/plink2_18_April_2015/plink --bfile ${outpath}/plink_cbr2_chr${CHR}_recode_pruned --extract plink_cbr2_3_4_common_variant.list --make-bed --out ${outpath}/plink_cbr2_chr${CHR}_recode_pruned_common
/nfs/team151/software/plink2_18_April_2015/plink --bfile ${outpath}/plink_cbr2_chr${CHR}_recode_pruned_common --freq --out ${outpath}/plink_cbr2_chr${CHR}_recode_pruned_freq

cbr3_imp_path=/lustre/scratch114/teams/soranzo/projects/Clean_Imputed_data/cbr_phase3/imputed/uk10k1kg.shapeit
/nfs/team151/software/plink2_18_April_2015/plink --gen ${cbr3_imp_path}/chr${CHR}.gen.gz --sample ${cbr3_imp_path}/chr${CHR}.01.gen_samples --oxford-single-chr ${CHR} --indep-pairwise 5000 1000 0.1 --out ${outpath}/plink_cbr3_chr${CHR}_recode
/nfs/team151/software/plink2_18_April_2015/plink --gen ${cbr3_imp_path}/chr${CHR}.gen.gz --sample ${cbr3_imp_path}/chr${CHR}.01.gen_samples --oxford-single-chr ${CHR} --extract plink_cbr3_chr${CHR}_recode.prune.in --make-bed --out ${outpath}/plink_cbr3_chr${CHR}_recode_pruned
/nfs/team151/software/plink2_18_April_2015/plink --bfile ${outpath}/plink_cbr3_chr${CHR}_recode_pruned --extract plink_cbr2_3_4_common_variant.list --make-bed --out ${outpath}/plink_cbr3_chr${CHR}_recode_pruned_common
/nfs/team151/software/plink2_18_April_2015/plink --bfile ${outpath}/plink_cbr3_chr${CHR}_recode_pruned_common --freq --out ${outpath}/plink_cbr3_chr${CHR}_recode_pruned_freq
# /nfs/team151/software/plink2_18_April_2015/plink --gen ${cbr3_imp_path}/chr1.gen.gz --sample ${cbr3_imp_path}/chr1.01.gen_samples --oxford-single-chr 1 --extract plink_cbr3_chr1_recode.prune.in --make-bed --out plink_cbr3_chr1_recode_pruned

cbr4_imp_path=/lustre/scratch114/teams/soranzo/projects/Clean_Imputed_data/cbr_phase4/imputed/uk10k1kg.shapeit
/nfs/team151/software/plink2_18_April_2015/plink --gen ${cbr4_imp_path}/chr${CHR}.gen.gz --sample ${cbr4_imp_path}/chr${CHR}.01.gen_samples --oxford-single-chr ${CHR} --indep-pairwise 5000 1000 0.1 --out ${outpath}/plink_cbr4_chr${CHR}_recode
/nfs/team151/software/plink2_18_April_2015/plink --gen ${cbr4_imp_path}/chr${CHR}.gen.gz --sample ${cbr4_imp_path}/chr${CHR}.01.gen_samples --oxford-single-chr ${CHR} --extract plink_cbr4_chr${CHR}_recode.prune.in --make-bed --out ${outpath}/plink_cbr4_chr${CHR}_recode_pruned
/nfs/team151/software/plink2_18_April_2015/plink --bfile ${outpath}/plink_cbr4_chr${CHR}_recode_pruned --extract plink_cbr2_3_4_common_variant.list --make-bed --out ${outpath}/plink_cbr4_chr${CHR}_recode_pruned_common
/nfs/team151/software/plink2_18_April_2015/plink --bfile ${outpath}/plink_cbr4_chr${CHR}_recode_pruned_common --freq --out ${outpath}/plink_cbr4_chr${CHR}_recode_pruned_freq

#fix missing id problem
awk '{if($2==".") print $1,"chr"$1":"$4,$3,$4,$5,$6;else print $0}' ${outpath}/plink_cbr2_chr${CHR}_recode_pruned_common.bim | tr " " "\t" > ${outpath}/plink_cbr2_chr${CHR}_recode_pruned_common_mod.bim
mv ${outpath}/plink_cbr2_chr${CHR}_recode_pruned_common.bim ${outpath}/plink_cbr2_chr${CHR}_recode_pruned_common_OLD.bim
mv ${outpath}/plink_cbr2_chr${CHR}_recode_pruned_common_mod.bim ${outpath}/plink_cbr2_chr${CHR}_recode_pruned_common.bim
awk '{if($2==".") print $1,"chr"$1":"$4,$3,$4,$5,$6;else print $0}' ${outpath}/plink_cbr3_chr${CHR}_recode_pruned_common.bim | tr " " "\t" > ${outpath}/plink_cbr3_chr${CHR}_recode_pruned_common_mod.bim
mv ${outpath}/plink_cbr3_chr${CHR}_recode_pruned_common.bim ${outpath}/plink_cbr3_chr${CHR}_recode_pruned_common_OLD.bim
mv ${outpath}/plink_cbr3_chr${CHR}_recode_pruned_common_mod.bim ${outpath}/plink_cbr3_chr${CHR}_recode_pruned_common.bim
awk '{if($2==".") print $1,"chr"$1":"$4,$3,$4,$5,$6;else print $0}' ${outpath}/plink_cbr4_chr${CHR}_recode_pruned_common.bim | tr " " "\t" > ${outpath}/plink_cbr4_chr${CHR}_recode_pruned_common_mod.bim
mv ${outpath}/plink_cbr4_chr${CHR}_recode_pruned_common.bim ${outpath}/plink_cbr4_chr${CHR}_recode_pruned_common_OLD.bim
mv ${outpath}/plink_cbr4_chr${CHR}_recode_pruned_common_mod.bim ${outpath}/plink_cbr4_chr${CHR}_recode_pruned_common.bim



#merge data together
# first remove allegedly multiallelic sites
/nfs/team151/software/plink2_18_April_2015/plink --bfile ${outpath}/plink_cbr2_chr${CHR}_recode_pruned_common --exclude plink_cbrALL_chr1_recode_pruned_MULTIALL.missnp --make-bed --out ${outpath}/plink_cbr2_chr${CHR}_recode_pruned_cleaned
/nfs/team151/software/plink2_18_April_2015/plink --bfile ${outpath}/plink_cbr3_chr${CHR}_recode_pruned_common --exclude plink_cbrALL_chr1_recode_pruned_MULTIALL.missnp --make-bed --out ${outpath}/plink_cbr3_chr${CHR}_recode_pruned_cleaned
/nfs/team151/software/plink2_18_April_2015/plink --bfile ${outpath}/plink_cbr4_chr${CHR}_recode_pruned_common --exclude plink_cbrALL_chr1_recode_pruned_MULTIALL.missnp --make-bed --out ${outpath}/plink_cbr4_chr${CHR}_recode_pruned_cleaned



/nfs/team151/software/plink2_18_April_2015/plink --bfile ${outpath}/plink_cbr2_chr${CHR}_recode_pruned_cleaned --merge-list ${outpath}/merge_list.list --make-bed --out ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned 

# fgrep Warning ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned.log | cut -f 3,5 -d " "| sed "s/'//g"|tr " " "\n" > ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned_duppos.list

/nfs/team151/software/plink2_18_April_2015/plink --bfile ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned  --exclude ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned_duppos.list --geno --mind --make-bed --out ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned_nodup
/nfs/team151/software/plink2_18_April_2015/plink --bfile ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned  --remove /lustre/scratch113/teams/soranzo/users/mc14/HE2checks/cbr4_overlaps_fam.txt --geno --mind --make-bed --out ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned_nodup_nodupsamples


#flash pca
# /nfs/team151/software/flashpca/flashpca --bfile plink_cbr3_chr1_recode_pruned --suffix _plink_cbr3_chr1_recode_pruned.txt

# /nfs/team151/software/flashpca/flashpca --bfile ${outpath}/plink_cbr2_chr${CHR}_recode_pruned --suffix _plink_cbr2_chr${CHR}_recode_pruned.txt
# /nfs/team151/software/flashpca/flashpca --bfile ${outpath}/plink_cbr3_chr${CHR}_recode_pruned --suffix _plink_cbr3_chr${CHR}_recode_pruned.txt
# /nfs/team151/software/flashpca/flashpca --bfile ${outpath}/plink_cbr4_chr${CHR}_recode_pruned --suffix _plink_cbr4_chr${CHR}_recode_pruned.txt

# pca with KING

echo "/nfs/team151/software/bin/king -b ${outpath}/plink_cbr2_chr${CHR}_recode_pruned.bed --mds --ibs --prefix plink_cbr2_chr${CHR}_recode_pruned_mds" | bsub -J"kinship" -o"%J_kinship.o" -M 5000 -R"select[mem>=5000] rusage[mem=5000]" -q yesterday     
echo "/nfs/team151/software/bin/king -b ${outpath}/plink_cbr3_chr${CHR}_recode_pruned.bed --mds --ibs --prefix plink_cbr3_chr${CHR}_recode_pruned_mds" | bsub -J"kinship" -o"%J_kinship.o" -M 5000 -R"select[mem>=5000] rusage[mem=5000]" -q yesterday     
echo "/nfs/team151/software/bin/king -b ${outpath}/plink_cbr4_chr${CHR}_recode_pruned.bed --mds --ibs --prefix plink_cbr4_chr${CHR}_recode_pruned_mds" | bsub -J"kinship" -o"%J_kinship.o" -M 5000 -R"select[mem>=5000] rusage[mem=5000]" -q yesterday     

echo "/nfs/team151/software/bin/king -b ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned.bed --mds --ibs --prefix ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned_mds" | bsub -J"kinship" -o"%J_kinship.o" -M 5000 -R"select[mem>=5000] rusage[mem=5000]" -q yesterday     
echo "/nfs/team151/software/bin/king -b ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned_nodup.bed --mds --ibs --prefix ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned_nodup_mds" | bsub -J"kinship" -o"%J_kinship.o" -M 5000 -R"select[mem>=5000] rusage[mem=5000]" -q yesterday     

echo "/nfs/team151/software/bin/king -b ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned.bed --pca 10 --prefix ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned_pca" | bsub -J"kinship" -o"%J_kinship.o" -M 5000 -R"select[mem>=5000] rusage[mem=5000]" -q yesterday     
# plink_cbr2_chr1_recode_pruned_mdspc.ped
echo "/nfs/team151/software/bin/king -b ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned_nodup_nodupsamples.bed --mds --ibs --prefix ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned_nodup_nodupsamples_mds" | bsub -J"kinship" -o"%J_kinship.o" -M 5000 -R"select[mem>=5000] rusage[mem=5000]" -q yesterday     

echo "/nfs/team151/software/bin/king -b ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned_nodup_nodupsamples.bed --mds --prefix ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned_nodup_nodupsamples_mds_noibs" | bsub -J"kinship" -o"%J_kinship.o" -M 5000 -R"select[mem>=5000] rusage[mem=5000]" -q yesterday     
# plink_cbr3_chr1_recode_pruned_mdspc.ped
# plink_cbr4_chr1_recode_pruned_mdspc.ped
echo "/nfs/team151/software/bin/king -b ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned_nodup_nodupsamples.bed --pca 10 --prefix ${outpath}/plink_cbrALL_chr${CHR}_recode_pruned_cleaned_nodupsamples_pca" | bsub -J"kinship" -o"%J_kinship.o" -M 5000 -R"select[mem>=5000] rusage[mem=5000]" -q yesterday     

# #PLOT PCA with R
# # cbr2_pca <- read.table("plink_cbr2_chr1_recode_pruned_mdspc.ped", header=F)
# # cbr3_pca <- read.table("plink_cbr3_chr1_recode_pruned_mdspc.ped", header=F)
# # cbr4_pca <- read.table("plink_cbr4_chr1_recode_pruned_mdspc.ped", header=F)

# plink_cbrALL_chr1_recode_pruned_cleaned_mdspc.ped
# plink_cbrALL_chr1_recode_pruned_cleaned_mdspc.dat

# 7159267_kinship.o
# plink_cbrALL_chr1_recode_pruned_cleaned_nodup_nodupsamples_mdspc.ped
# plink_cbrALL_chr1_recode_pruned_cleaned_nodup_nodupsamples_mdspc.dat

# 7165970_kinship.o
# plink_cbrALL_chr1_recode_pruned_cleaned_pcapc.ped
# plink_cbrALL_chr1_recode_pruned_cleaned_pcapc.dat

cbrALL_pca <- read.table("plink_cbrALL_chr1_recode_pruned_cleaned_nodup_mdspc.ped", header=F)
# cbrALL_pca <- read.table("plink_cbrALL_chr1_recode_pruned_cleaned_nodup_nodupsamples_mdspc.ped", header=F)

cbrALL_pca <- read.table("plink_cbrALL_chr1_recode_pruned_cleaned_pcapc.ped", header=F)
# cbrALL_pca <- read.table("plink_cbrALL_chr1_recode_pruned_cleaned_mdspc.ped", header=F)
cbrALL_pca$V2 <- NULL
cbrALL_pca$V3 <- NULL
cbrALL_pca$V4 <- NULL
cbrALL_pca$V5 <- NULL
cbrALL_pca$V6 <- NULL
# cbr2_pca$V2 <- NULL
# cbr3_pca$V2 <- NULL
# cbr4_pca$V2 <- NULL
# cbr2_pca$V3 <- NULL
# cbr3_pca$V3 <- NULL
# cbr4_pca$V3 <- NULL
# cbr2_pca$V4 <- NULL
# cbr3_pca$V4 <- NULL
# cbr4_pca$V4 <- NULL
# cbr2_pca$V5 <- NULL
# cbr3_pca$V5 <- NULL
# cbr4_pca$V5 <- NULL
# cbr2_pca$V6 <- NULL
# cbr3_pca$V6 <- NULL
# cbr4_pca$V6 <- NULL

# colnames(cbr2_pca) <- c("ID",comp)
# colnames(cbr3_pca) <- c("ID",comp)
# colnames(cbr4_pca) <- c("ID",comp)
comp <- paste(rep("C",times=20),seq(1,20,1),sep="")
colnames(cbrALL_pca) <- c("ID",comp)

#fam files for all sets
cbr_fam2 <- read.table("plink_cbr2_chr1_recode_pruned_common.fam",header=F)
cbr_fam3 <- read.table("plink_cbr3_chr1_recode_pruned_common.fam",header=F)
cbr_fam4 <- read.table("plink_cbr4_chr1_recode_pruned_common.fam",header=F)

cbr_fam2$cohort <- "cbr2"
cbr_fam3$cohort <- "cbr3"
cbr_fam4$cohort <- "cbr4"

cbr_famALL <- rbind(cbr_fam2,cbr_fam3,cbr_fam4)
cbr_famALL$V2 <- NULL
cbr_famALL$V3 <- NULL
cbr_famALL$V4 <- NULL
cbr_famALL$V5 <- NULL
cbr_famALL$V6 <- NULL

colnames(cbr_famALL) <- c("ID","cohort")

cbrALL_pca_cohort <- merge(cbrALL_pca,cbr_famALL,by="ID", all=TRUE)

# ###########################################################################
# # consecutive components
# cohorts <- c("cbr2","cbr3","cbr4")
cohorts <- c("cbrALL")
for (cohort in cohorts){
# cohort <- cohorts[1]
	
current_set <- get(paste(cohort,"_pca_cohort",sep=""))
# current_set <- get(paste(cohort,"_pca",sep=""))
xmax <- max(current_set$C1)-0.001
ymax <- max(current_set$C2)+0.001
basefolder <- getwd()
# pdf(paste(basefolder,"PCA_",cohort,"_batch_PRUNED_",type,"_",filter,"_gwa_cons.pdf",sep=""),width=11.7, height=8.3)
# jpeg(paste(basefolder,"/PCA_",cohort,"_cons_nodupsamples_mdspc.jpg",sep=""),width=1754, height=1024,pointsize = 20)
jpeg(paste(basefolder,"/PCA_",cohort,"_cons_pca_col.jpg",sep=""),width=1754, height=1024,pointsize = 20)
par(lab=c(4,4,6),mfrow=c(2,3))
plot(current_set[which(current_set$cohort == "cbr2"),]$C1,current_set[which(current_set$cohort == "cbr2"),]$C2,xlim=c(min(current_set$C1)-0.001,max(current_set$C1)+0.001),ylim=c(min(current_set$C2)-0.001,max(current_set$C2)+0.001),main=paste("PCA C1 vs C2, ",cohort,sep=""),xlab="C1",ylab="C2")
points(current_set[which(current_set$cohort == "cbr3"),]$C1,current_set[which(current_set$cohort == "cbr3"),]$C2,col="red")
points(current_set[which(current_set$cohort == "cbr4"),]$C1,current_set[which(current_set$cohort == "cbr4"),]$C2,col="green")
# legend(xmax,ymax,legend=c("BGI","SC"),fill=c("black","red"))

plot(current_set[which(current_set$cohort == "cbr2"),]$C3,current_set[which(current_set$cohort == "cbr2"),]$C4,xlim=c(min(current_set$C3),max(current_set$C3)),ylim=c(min(current_set$C4),max(current_set$C4)),main=paste("PCA C3 vs C4, ",cohort,sep=""),xlab="C3",ylab="C4")
# points(current_set[which(current_set$SEQ == 2),]$C3,current_set[which(current_set$SEQ == 2),]$C4,col="red")
points(current_set[which(current_set$cohort == "cbr3"),]$C3,current_set[which(current_set$cohort == "cbr3"),]$C4,col="red")
points(current_set[which(current_set$cohort == "cbr4"),]$C3,current_set[which(current_set$cohort == "cbr4"),]$C4,col="green")

plot(current_set[which(current_set$cohort == "cbr2"),]$C5,current_set[which(current_set$cohort == "cbr2"),]$C6,xlim=c(min(current_set$C5),max(current_set$C5)),ylim=c(min(current_set$C6),max(current_set$C6)),main=paste("PCA C5 vs C6, ",cohort,sep=""),xlab="C5",ylab="C6")
# points(current_set[which(current_set$SEQ == 2),]$C5,current_set[which(current_set$SEQ == 2),]$C6,col="red")
points(current_set[which(current_set$cohort == "cbr3"),]$C5,current_set[which(current_set$cohort == "cbr3"),]$C6,col="red")
points(current_set[which(current_set$cohort == "cbr4"),]$C5,current_set[which(current_set$cohort == "cbr4"),]$C6,col="green")

plot(current_set[which(current_set$cohort == "cbr2"),]$C7,current_set[which(current_set$cohort == "cbr2"),]$C8,xlim=c(min(current_set$C7),max(current_set$C7)),ylim=c(min(current_set$C8),max(current_set$C8)),main=paste("PCA C7 vs C8, ",cohort,sep=""),xlab="C7",ylab="C8")
# points(current_set[which(current_set$SEQ == 2),]$C7,current_set[which(current_set$SEQ == 2),]$C8,col="red")
points(current_set[which(current_set$cohort == "cbr3"),]$C7,current_set[which(current_set$cohort == "cbr3"),]$C8,col="red")
points(current_set[which(current_set$cohort == "cbr4"),]$C7,current_set[which(current_set$cohort == "cbr4"),]$C8,col="green")

plot(current_set[which(current_set$cohort == "cbr2"),]$C9,current_set[which(current_set$cohort == "cbr2"),]$C10,xlim=c(min(current_set$C9),max(current_set$C9)),ylim=c(min(current_set$C10),max(current_set$C10)),main=paste("PCA C9 vs C10, ",cohort,sep=""),xlab="C9",ylab="C10")
# points(current_set[which(current_set$SEQ == 2),]$C9,current_set[which(current_set$SEQ == 2),]$C10,col="red")
points(current_set[which(current_set$cohort == "cbr3"),]$C9,current_set[which(current_set$cohort == "cbr3"),]$C10,col="red")
points(current_set[which(current_set$cohort == "cbr4"),]$C9,current_set[which(current_set$cohort == "cbr4"),]$C10,col="green")

dev.off()

}
# # plot(cbr2_pca$C1,cbr2_pca$C2,col="red")
# # plot(cbr3_pca$C1,cbr3_pca$C2,col="green")
# # plot(cbr4_pca$C1,cbr4_pca$C2,col="blue")

#fix gender issues
#temporary id samples
paste -d " " <(tail -n+3 chr5.32.gen_samples) <(tail -n+3 chr5.32.gen_samples| cut -f 1 -d " "| cut -f 3 -d "_") > cbr4.gen_samples.tmp

awk 'FNR==NR { a[$8]=$0; next } $3 in a { print $0,a[$8] }' /lustre/scratch114/teams/soranzo/projects/Clean_Imputed_data/cbr_phase4/imputed/uk10k1kg.shapeit/cbr4.gen_samples.tmp /lustre/scratch113/teams/soranzo/users/mc14/HE2checks/cbr4_samples_sex.list