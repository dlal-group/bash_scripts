#!/usr/local/bin/bash

#Pipeline for batch effect analysis check
#awk '{print $1,$2}' data.fam | split -d -a 3 -l 100 - tmp.list #split the list of samples
cohort=$1
indir=$2
genopath=$3
freqpath=$4
outdir=$5

# maf=0.01 #only common variants
mkdir -p ${outdir}

a=$[`ls ${indir}/tmp.list*| wc -l`-1]

let i=0;let j=0;

while [ $i -le $a ]
do
    while [ $j -le $a ]
    do

        bsub -J "plink_run_${i}_${j}" -o "${outdir}/%J_plink_run_${i}_${j}.o" -M5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- plink --noweb \
        --bfile ${genopath} \
        --remove ~/fvg_seq/REL-2014-01-09/v1_stats/PLINK/GWAS_OVERLAP/single_sample_check.list \
        --read-freq ${freqpath} \
        --genome \
        --genome-lists ${indir}/tmp.list`printf "%03i\n" $i` \
        ${indir}/tmp.list`printf "%03i\n" $j` \
        --out ${outdir}/data.sub.$i.$j

        let j=$j+1
    done

    let i=$i+1
    let j=$i

done

echo "cat ${outdir}/data.sub*genome > ${outdir}/data.genome;rm ${outdir}/data.sub.*" | bsub -J "concat_mds" -o "${outdir}/%J_concat_mds.o" -M2000 -w "ended(plink_run_*)" -R"select[mem>=2000] rusage[mem=2000]"

# ##4. calculate MDS

# maf=0.01 #only common variants
# #max_maf=0.05
mkdir -p ${outdir}/MDS

bsub -J "plink_mds_all" -o "${outdir}/%J_plink_mds_all.o" -M8000 -w "ended(concat_mds)" -R"select[mem>=8000] rusage[mem=8000]" -q normal -- plink --noweb \
--bfile ${genopath} \
--read-genome ${outdir}/data.genome \
--cluster --mds-plot 10 \
--out ${outdir}/MDS/all_chr_merged

#####################################
###R code to plot with FVG data


rm(list=ls())
cohort <- "FVG"
fvg <- read.table("/nfs/users/nfs_m/mc14/Work/SANGER/FVG/SEQ_CALLING/DOCS/samples_resume.txt",header=T)
mds <- read.table("all_chr_merged.mds",header=T)
#mds <- read.table("chr20.mds.mod",sep="\t",header=T)
merged_new <- merge(mds,fvg,by.x="IID",by.y="clinic_id",all.x)
villages <- unique(merged_new$village)
#remove the sample with all missing,otherwise the plot for c1/c2 is too shrinked
merged_new <- merged_new[-which(merged_new$IID == "591350"),]
#plot different villages
jpeg(paste("PCA_",cohort,".jpg",sep=""),width=1754, height=1024,pointsize = 20)
par(lab=c(4,4,6),mfrow=c(2,3))
for (c in seq(4,which(colnames(merged_new) == "C10"),by=2)) {
    xmax <- max(merged_new[c])
    ymax <- max(merged_new[c+1])
    plot(merged_new[which(merged_new$village == villages[1]),c],merged_new[which(merged_new$village == villages[1]),c+1],xlim=c(min(merged_new[c])-0.001,max(merged_new[c])+0.001),ylim=c(min(merged_new[c+1])-0.001,max(merged_new[c+1])+0.001),main=paste("PCA ",colnames(merged_new)[c]," vs ",colnames(merged_new)[c+1]," ",cohort,sep=""),xlab=colnames(merged_new)[c],ylab=colnames(merged_new)[c+1], col=colors()[72])
    vill_cols <- colors()[72]
    for ( i in 2:length(villages)){
        points(merged_new[which(merged_new$village == villages[i]),c],merged_new[which(merged_new$village == villages[i]),c+1],col=colors()[72+(i*5)])
        vill_cols <- c(vill_cols,colors()[72+(i*5)])
    }
}
plot.new()
legend("center",legend=villages,fill=vill_cols,xpd=TRUE,cex=2, bty="n", title="Villages:")

dev.off()

jpeg(paste("PCA_cons_",cohort,".jpg",sep=""),width=1754, height=1024,pointsize = 20)
par(lab=c(4,4,6),mfrow=c(2,3))
for (c in seq(4,which(colnames(merged_new) == "C4"),by=1)) {
    xmax <- max(merged_new[c])
    ymax <- max(merged_new[c+1])
    plot(merged_new[which(merged_new$village == villages[1]),c],merged_new[which(merged_new$village == villages[1]),c+1],xlim=c(min(merged_new[c])-0.001,max(merged_new[c])+0.001),ylim=c(min(merged_new[c+1])-0.001,max(merged_new[c+1])+0.001),main=paste("PCA ",colnames(merged_new)[c]," vs ",colnames(merged_new)[c+1]," ",cohort,sep=""),xlab=colnames(merged_new)[c],ylab=colnames(merged_new)[c+1], col=colors()[72])
    vill_cols <- colors()[72]
    for ( i in 2:length(villages)){
        points(merged_new[which(merged_new$village == villages[i]),c],merged_new[which(merged_new$village == villages[i]),c+1],col=colors()[72+(i*5)])
        vill_cols <- c(vill_cols,colors()[72+(i*5)])
    }
}
plot.new()
legend("center",legend=villages,fill=vill_cols,xpd=TRUE,cex=2, bty="n", title="Villages:")

dev.off()
