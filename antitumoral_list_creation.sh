#!/usr/local/bin/bash
#Script used to create merged files for the population
POP=$1 #FVG
markerlist=$2 #/nfs/users/nfs_m/mc14/Work/SANGER/FVG/ANTI_TUMORAL_DRUGS/BREAST_REGIONS/FREQ/SNPs/antitumoral_marker_pop_freq.csv
outpath=$3 #/nfs/users/nfs_m/mc14/Work/SANGER/FVG/ANTI_TUMORAL_DRUGS/BREAST_REGIONS/FREQ/SNPs
outfile=$4 #ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.antitumoral_marker.vcf

(zcat /lustre/scratch113/projects/fvg_seq/variant_refinemet/annotations/1TGP/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz | head -1000| grep "^#";while read -a line
do

   chr=${line[0]}
   site=${line[1]}
   tabix /lustre/scratch113/projects/fvg_seq/variant_refinemet/annotations/1TGP/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz ${chr}| fgrep -w ${site}

done < <(uniq ${markerlist}) | sort -g -k1,1 -k2,2) | bgzip -c > ${outpath}/${outfile}.gz

tabix -p vcf ${outpath}/${outfile}.gz
#
##extract freq data in tab format
(echo "CHROM POS ID REF ALT AN AC AF AMR_AF ASN_AF AFR_AF EUR_AF";bcftools2 query ${outpath}/${outfile}.gz -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\t%INFO/AF\t%INFO/AMR_AF\t%INFO/ASN_AF\t%INFO/AFR_AF\t%INFO/EUR_AF\n") | tr "\t" " " > ${outpath}/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.antitumoral_marker.csv

#merged pop data with 1000G data
#Case 1 : plink freq files
(echo "CHROM POS ID REF ALT AN AC AF AMR_AF ASN_AF AFR_AF EUR_AF ${POP}_A1 ${POP}_A2 ${POP}_MAF ${POP}_AF ${POP}_NCHROBS";awk 'FNR==NR{a[$3]=$0;next}{ print a[$2],$3,$4,$5,1-$5,$6}' <(tail -n+2 ${outpath}/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.antitumoral_marker.csv) <(tail -n+2 ${markerlist}))|tr "\t" " " > ${outpath}/antitumoral_marker_pop_freq_merged.csv

#add columns with check allele
(echo "CHROM POS ID REF ALT AN AC AF AMR_AF ASN_AF AFR_AF EUR_AF ${POP}_A1 ${POP}_A2 ${POP}_MAF ${POP}_AF ${POP}_NCHROBS ${POP}_CHK_ALL ${POP}_CHK_AF";awk '{if($5==$13) print $0,$13,$15;else print $0,$14,$16}' <(tail -n+2 ${outpath}/antitumoral_marker_pop_freq_merged.csv)) > ${outpath}/antitumoral_marker_pop_freq_merged_complete.csv
#Case 2 : imputation result
# (echo "CHROM POS ID REF ALT AN AC AF AMR_AF ASN_AF AFR_AF EUR_AF ${POP}_A1 ${POP}_A2 ${POP}_MAF ${POP}_AF ${POP}_NCHROBS";awk 'FNR==NR{a[$3]=$0;next}{ print a[$2],$3,$4,$5,1-$5,$6}' <(tail -n+2 ${outpath}/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.antitumoral_marker.csv) <(tail -n+2 ${markerlist}))|tr "\t" " " > ${outpath}/antitumoral_marker_pop_freq_merged.csv


# (echo "CHROM POS ID REF ALT AN AC AF AMR_AF ASN_AF AFR_AF EUR_AF FVG_A1 FVG_A2 FVG_MAF FVG_AF FVG_NCHROBS FVG_CHK_ALL FVG_CHK_AF";awk '{if($5==$13) print $0,$13,$15;else print $0,$14,$16}' <(tail -n+2 antitumoral_marker_pop_freq_merged.csv)) > antitumoral_marker_pop_freq_merged_complete.csv