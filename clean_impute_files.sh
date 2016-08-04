#!/usr/bin/env bash
#

pop=$1
chr=$2
basefolder=$3
# basefolder=/netapp02/data/imputation/INGI_TGP3/impute
# for pop in CARL FVG VBI
# do

# for chr in {1..22}
# do
# pop=FVG
# chr=1
last_chunk=`ls ${basefolder}/${pop}/${chr}/*.gz| sort| tail -n1`
last_chunk_n=`echo ${last_chunk}| cut -f 2 -d "."`
#first, merge all chunks for gen and info
for chunk in `seq $last_chunk_n`
do
chunkStr=`printf "%02d" $chunk`
echo ${chr} ${chunkStr}

cat ${basefolder}/${pop}/${chr}/chr${chr}.${chunkStr}.gen.gz >> ${basefolder}/${pop}/chr${chr}.gen_tmp1.gz
cat ${basefolder}/${pop}/${chr}/chr${chr}.${chunkStr}.gen_info | sed 's,'"/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/CARL_FVG_VBI_TGP3_ALL/${chr}/${chr}.INGI_REF.CARL_FVG_VBI_TGP3_ALL.*.legend.gz:"',,g' >> ${basefolder}/${pop}/chr${chr}.gen_tmp1_info

done

#create the COMPLETE merged imputed files
echo "Creating merged complete files...."
mkdir -p ${basefolder}/${pop}/MERGED/ALL
mkdir -p ${basefolder}/${pop}/MERGED/CLEANED

zcat ${basefolder}/${pop}/chr${chr}.gen_tmp1.gz | sed 's,'"/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/CARL_FVG_VBI_TGP3_ALL/${chr}/${chr}.INGI_REF.CARL_FVG_VBI_TGP3_ALL.*.legend.gz:"',,g' > ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen.gz
(echo "snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0";fgrep -v -w rs_id ${basefolder}/${pop}/chr${chr}.gen_tmp1_info) > ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen_info

echo "Cleaning step...."
#than clean up useless text and recode long alleles, creating a map file for decoding
zcat ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen.gz | awk '{if (x[$3]) { x_count[$3]++; print $0; if (x_count[$3] == 0) { print x[$3] } } x[$3] = $0}' | gzip -c > ${basefolder}/${pop}/chr${chr}.gen_tmp2.gz
awk '{if (x[$3]) { x_count[$3]++; print $0; if (x_count[$3] == 0) { print x[$3] } } x[$3] = $0}' ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen_info  > ${basefolder}/${pop}/chr${chr}.gen_tmp2_info

#extract duplicates for imputation results:
#from info and geno files for pos and alleles length (cols 3,4,5)
# awk '{if (x[$3,length($4),length($5)]) { x_count[$3,length($4),length($5)]++; print $0; if (x_count[$3,length($4),length($5)] > 1) { print x[$3,length($4),length($5)] } } x[$3,length($4),length($5)] = $0}' 

# done
# done