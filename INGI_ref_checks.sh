#!/usr/local/bin/bash

#Script to prepare data for plotting after imputation
    # for pan in uk10k1kg.ref
# for pan in INGI.shapeit CARL.shapeit FVG.shapeit VBI.shapeit 1000Gph1.shapeit 1000GP_Phase3.shapeit INGI_1000GPh3.shapeit uk10k1kg.ref
    # for pan in 1000GP_Phase3.shapeit

# pop=CARL
    # for pan in INGI.shapeit CARL.shapeit FVG.shapeit VBI.shapeit 1000Gph1.shapeit INGI_1000GPh3.shapeit uk10k1kg.ref
    # pan=1000Gph1.shapeit
    # pan=uk10k1kg.ref
    # do
basefolder2="/lustre/scratch113/projects/esgi-vbseq/31032016_IMPUTATION"
basefolder="/lustre/scratch113/teams/soranzo/users/jh21/imputed"

# for pop in CARL FVG INCIPE2 VBI
for pop in carl fvg incipe2 vbi
# for pop in INCIPE2
do
# for pan in CARL_FVG_VBI.shapeit CARL_FVG_VBI_TSI.shapeit CARL.shapeit FVG.shapeit VBI.shapeit
for pan in uk10k1kg.ref
do

    mkdir -p ${basefolder2}/${pop^^}/${pan}
        for chr in 21
        do
            maf_bins="0,0.005,0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50"
            # chr=1
            (fgrep -h -v position ${basefolder}/${pop}/${pan}/chr${chr}.*.gen_info) | gzip -c > ${basefolder2}/${pop}/${pan}/chr${chr}.gen_info.gz
            # gzip ${pop}/${pan}/chr${chr}.gen_info
            (echo "CHROM RS_ID POS EXP_FREQ_A1 INFO TYPE INFO_TYPE0 CONCORD_TYPE0 r2_TYPE0 COHORT PANEL MAF BIN";(zgrep -v position ${basefolder2}/${pop}/${pan}/chr${chr}.gen_info.gz | cut -f 2,3,6,7,9- -d " "| awk -v chrom=$chr -v cohort=$pop -v panel=$pan '{if($3 <= 0.5 ) print chrom,$0,toupper(cohort),panel,$3; else print chrom,$0,toupper(cohort),panel,1-$3}'|awk -v bins=$maf_bins '
            {n=split(bins,mafs,",");}{
                for (i=1;i<=n;i++){
                    if (mafs[i]==0){
                        if($(NF) <= mafs[i]){
                            print $0,mafs[i]
                        }
                    }else{
                        if($(NF) <= mafs[i] && $(NF) > mafs[i-1]){
                            print $0,mafs[i]
                        }

                    }
                }
            }'))| gzip -c > ${basefolder2}/${pop}/${pan}/chr${chr}.gen_info_partial_t2.gz

            # (zcat /lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/05272015_MERGED_REF_PANEL/IMPUTED/${pop}/${pan}/chr${chr}.gen_info_partial.gz| head -1;(zcat /lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/05272015_MERGED_REF_PANEL/IMPUTED/${pop}/${pan}/chr${chr}.gen_info_partial.gz| awk '$6 == 2')) | gzip -c > /lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/05272015_MERGED_REF_PANEL/IMPUTED/${pop}/${pan}/chr${chr}.gen_info_partial_t2.gz 
        done
    done
done
