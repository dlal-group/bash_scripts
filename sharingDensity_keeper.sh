#!/usr/local/bin/bash

reg_file=$1 #format:${MATCH}.${minDens}.shareDens_R*.to_include

MAP=$2 #PATH of plink format files in cM
# MAP=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/

#modify to extract start and end snp to extract region from plink
# for reg_file in `ls ${MATCH}.${minDens}.shareDens_R*.to_include`
# do
    # start_r=`head -1 ${reg_file} | awk '{print $1*0.5}'`
    # end_r=`tail -1 ${reg_file} | awk '{print $1*0.5}'`
    start_r=`head -1 ${reg_file} | cut -f 1`
    end_r=`tail -1 ${reg_file} | cut -f 1`
    CHR=`echo ${reg_file} | cut -f 2 -d "."`
    reg_start=`awk -v rstrt=${start_r} -v rnd=$end_r '{if($3 >= rstrt && $3<= rnd ) print $2}' $MAP/ALL.${CHR}.non_missing.map |head -1`
    reg_end=`awk -v rstrt=${start_r} -v rnd=$end_r '{if($3 >= rstrt && $3<= rnd ) print $2}' $MAP/ALL.${CHR}.non_missing.map |tail -1`
    
    echo -e "Arguments and parameters:
    \rMatch file-> $MATCH
    \rMap file -> ALL.${CHR}.non_missing.map"

    echo "${reg_start} ${reg_end}" > ${reg_file}.keepsnps
# done

# fgrep -v -f $MATCH.excluded_regions $MATCH
# done < <(awk -v max_sd=${max_d} '$2>=max_sd' CEU.22.non_missing.match.shareDens) > excluded_regions.txt

# MATCH=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/IBD/GERMLINE/CHR22/CEU.22.non_missing.match
# MAP=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/CEU.22.non_missing.map
# start_r=`awk -v max_sd=${max_d} '$2>=max_sd' $MATCH.shareDens | awk '{print $1*0.5}'`
# end_r=`awk -v max_sd=${max_d} '$2>=max_sd' $MATCH.shareDens | awk '{print ($1+0.5)*0.5}'`

# avg=`awk '{sum+=$2} END { printf "%.8f", sum/NR}' CEU.22.non_missing.match.shareDens`
# std=`awk -v mu=${avg} '{sumsq += ($2 - mu)^2} END { printf "%.8f", sqrt((sumsq)/NR)}' CEU.22.non_missing.match.shareDens` 
# start_r=`awk -v max_sd=${max_d} '$2>=max_sd' CEU.22.non_missing.match.shareDens | awk '{print $1*0.5}'`
# end_r=`awk -v max_sd=${max_d} '$2>=max_sd' CEU.22.non_missing.match.shareDens | awk '{print ($1+0.5)*0.5}'`


# numero medio coppie di individui vs size in cM -> lo si fa tra unrelated e tra related e vedere cosa succede
