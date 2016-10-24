#!/usr/bin/env bash

# 24/10/2016
# Extract numbers for annotations with freqs and AC based on CADD thresholds:
# [0-5) 
# [5-15)
# [15-20)
# >= 20

# want to know how many for each threshold, what kind of consequences, and how many singleton
file=$1
v_type=$2
base_dir=`dirname ${file}`
file_name=`basename ${file}`

CADD=("${@:3}")

CADD_val=`echo ${#CADD[@]}`

mkdir -p ${base_dir}/CADD_STRAT

if [[ ${CADD_val} -eq 2 ]]; then
	#statements
	CADD_1=${CADD[0]}
	CADD_2=${CADD[1]}
	bcftools norm -m - ${file} | bcftools +fill-AN-AC | bcftools query -i'INFO/CADD_PHRED>=${CADD_1} && INFO/CADD_PHRED<${CADD_2} && TYPE="${v_type}"' -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CADD\t%INFO/CSQ\t%AC\t%AN\n" | awk 'BEGIN{OFS="\t"}{print $0, $(NF-1)/$(NF)}'| awk 'BEGIN{OFS="\t"}{if($(NF)<= 0.5) print $0, $(NF);else print $0, 1-$(NF)}' > ${base_dir}/CADD_STRAT/${file_name}.${v_type}.${CADD_1}_${CADD_2}.freq.tab
	
	#we need to split each table by LOF category
	for cat in frameshift splice_acceptor splice_donor stop_gain stop_loss
	do
		fgrep ${cat} ${base_dir}/CADD_STRAT/${file_name}.${vtype}.${CADD_1}_${CADD_2}.freq.tab > ${base_dir}/CADD_STRAT/${file_name}.${v_type}.${CADD_1}_${CADD_2}.${cat}.tab
	done

	bcftools norm -m - ${file} | bcftools +fill-AN-AC | bcftools view -G -i'INFO/CADD_PHRED>=${CADD_1} && INFO/CADD_PHRED<${CADD_2} && TYPE="${v_type}"' -O z -o ${base_dir}/CADD_STRAT/${file_name}.${v_type}.${CADD_1}_${CADD_2}.vcf.gz
	tabix -p vcf ${base_dir}/CADD_STRAT/${file_name}.${v_type}.${CADD_1}_${CADD_2}.vcf.gz
else
	CADD_1=${CADD[0]}
	bcftools norm -m - ${file} | bcftools +fill-AN-AC | bcftools query -i'INFO/CADD_PHRED>=${CADD_1} && TYPE="${v_type}"' -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CADD\t%INFO/CSQ\t%AC\t%AN\n" | awk 'BEGIN{OFS="\t"}{print $0, $(NF-1)/$(NF)}'| awk 'BEGIN{OFS="\t"}{if($(NF)<= 0.5) print $0, $(NF);else print $0, 1-$(NF)}' > ${base_dir}/CADD_STRAT/${file_name}.${v_type}.${CADD_1}.freq.tab
	
	#we need to split each table by LOF category
	for cat in frameshift splice_acceptor splice_donor stop_gain stop_loss
	do
		fgrep ${cat} ${base_dir}/CADD_STRAT/${file_name}.${v_type}.${CADD_1}.freq.tab > ${base_dir}/CADD_STRAT/${file_name}.${v_type}.${CADD_1}.${cat}.tab
	done

	bcftools norm -m - ${file} | bcftools +fill-AN-AC | bcftools view -G -i'INFO/CADD_PHRED>=${CADD_1} && TYPE="${v_type}"' -O z -o ${base_dir}/CADD_STRAT/${file_name}.${v_type}.${CADD_1}.vcf.gz
	tabix -p vcf ${base_dir}/CADD_STRAT/${file_name}.${v_type}.${CADD_1}.vcf.gz
fi
