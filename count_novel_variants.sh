#!/usr/bin/env bash

#chr is the prefix for the file with the whole path
mode=$1
 case $mode in
 	INGI )
		chr=$2
		bcftools norm -f /lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa -m - /lustre/scratch113/projects/esgi-vbseq/01032016_PANEL_SESOURCES/INGI/MERGE/${chr}.vcf.gz.norm.vcf.gz | bcftools +missing2ref -- -p | bcftools +fill-AN-AC | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%AN\n" | awk 'BEGIN{OFS="\t"}{print $0, $(NF-1)/$(NF)}'| awk 'BEGIN{OFS="\t"}{if($(NF)<= 0.5) print $0, $(NF);else print $0, 1-$(NF)}' >  /lustre/scratch113/projects/esgi-vbseq/01032016_PANEL_SESOURCES/INGI/MERGE/INGI.${chr}.MERGED.freq.tab
 		;;
	OUTBRED )
		chr=$2
		vtype=$3
		union_path=$4
		echo "${chr}"
		echo "${union_path}/${chr}/sites.txt.${vtype}.INGI.share_count.tab"
		echo "/lustre/scratch113/projects/esgi-vbseq/01032016_PANEL_SESOURCES/INGI/MERGE/${chr}.vcf.gz.norm.vcf.gz"
		echo "${union_path}/${chr}/INGI.${vtype}.novel.freq.tab"


		cat ${union_path}/${chr}/sites.txt.${vtype}.CARL_FVG.share_count.tab ${union_path}/${chr}/sites.txt.${vtype}.CARL.share_count.tab ${union_path}/${chr}/sites.txt.${vtype}.CARL_VBI_FVG.share_count.tab ${union_path}/${chr}/sites.txt.${vtype}.CARL_VBI.share_count.tab ${union_path}/${chr}/sites.txt.${vtype}.FVG.share_count.tab ${union_path}/${chr}/sites.txt.${vtype}.VBI_FVG.share_count.tab ${union_path}/${chr}/sites.txt.${vtype}.VBI.share_count.tab | sort -g -k2,2| uniq| tr " " "\t" > ${union_path}/${chr}/sites.txt.${vtype}.INGI.share_count.tab
		bcftools norm -f /lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa -m - -R ${union_path}/${chr}/sites.txt.${vtype}.INGI.share_count.tab /lustre/scratch113/projects/esgi-vbseq/01032016_PANEL_SESOURCES/INGI/MERGE/${chr}.vcf.gz.norm.vcf.gz | bcftools +missing2ref -- -p | bcftools +fill-AN-AC | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%AN\n" | awk 'BEGIN{OFS="\t"}{print $0, $(NF-1)/$(NF)}'| awk 'BEGIN{OFS="\t"}{if($(NF)<= 0.5) print $0, $(NF);else print $0, 1-$(NF)}' > ${union_path}/${chr}/INGI.${vtype}.novel.freq.tab
 		;;
 esac


