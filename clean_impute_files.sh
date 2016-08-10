#!/usr/bin/env bash
#

pop=$1
chr=$2
basefolder=$3
# basefolder=/netapp02/data/imputation/INGI_TGP3/impute
mode=$4

case $mode in
	STEP1 )
	echo "First step: MERGE and spurious text removal"
		last_chunk=`ls ${basefolder}/${pop}/${chr}/*.gz| sort| tail -n1`
		last_chunk_n=`echo ${last_chunk}| cut -f 2 -d "."`
		#first, merge all chunks for gen and info
		for chunk in `seq $last_chunk_n`
		do
		chunkStr=`printf "%02d" $chunk`
		echo ${chr} ${chunkStr}

		cat ${basefolder}/${pop}/${chr}/chr${chr}.${chunkStr}.gen.gz >> ${basefolder}/${pop}/chr${chr}.gen_tmp1.gz
		cat ${basefolder}/${pop}/${chr}/chr${chr}.${chunkStr}.gen_info | sed 's,'"/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/CARL_FVG_VBI_TGP3_ALL/${chr}/${chr}.INGI_REF.CARL_FVG_VBI_TGP3_ALL.*.legend.gz:"',,g' >> ${basefolder}/${pop}/chr${chr}.gen_tmp1_info

		#create the COMPLETE merged imputed files
		echo "Creating merged complete files...."
		mkdir -p ${basefolder}/${pop}/MERGED/ALL
		mkdir -p ${basefolder}/${pop}/MERGED/CLEANED
		mkdir -p ${basefolder}/${pop}/MERGED/CLEANED/RECODED

		zcat ${basefolder}/${pop}/chr${chr}.gen_tmp1.gz | sed 's,'"/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/CARL_FVG_VBI_TGP3_ALL/${chr}/${chr}.INGI_REF.CARL_FVG_VBI_TGP3_ALL.*.legend.gz:"',,g'| gzip -c > ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen.gz
		(echo "snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0";fgrep -v -w rs_id ${basefolder}/${pop}/chr${chr}.gen_tmp1_info) > ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen_info
		done
		;;
	STEP2 )
	echo "Second step: Multiallelic removal"
		#extract a list of sites to keep:this sites will be uniq by position!
		/home/cocca/scripts/bash_scripts/clean_impute_files.py ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen ${chr} ${pop} ${basefolder}/${pop}/MERGED/CLEANED

		#sort the keeplist
		sort -g -k2,2 ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen.to_keep -o ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen.to_keep
		#create the rsId keeplist for qctool
		cut -f 1 -d " " ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen.to_keep > ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen_rsID.to_keep
		#now we have the list of stuff we want to keep, so we'll just grep what we need
		(echo "snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0";fgrep -f ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen.to_keep ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen_info) > ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen_info
		qctool -g ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen.gz -incl-rsids ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen_rsID.to_keep -omit-chromosome -og ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen.gz 
		#create sample file
		cut -f -6 -d " " ${basefolder}/${pop}/${chr}/chr${chr}.05.gen_samples > ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen_samples
		;;
	STEP3 )
	echo "Third step: Long allele names recode"
		#and last thing we need is to RECODE long alleles names
		awk '{
		if (length($4)!=length($5)){
		if(length($4)>length($5)){
		printf $1" "$2" "$3" R D ";for(i=6;i<=NF;i++) printf "%s ", $i; printf "\n"
		} else {
		printf $1" "$2" "$3" R I ";for(i=6;i<=NF;i++) printf "%s ", $i; printf "\n"
		}
		} else {
		print $0
		}}' ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen_info > ${basefolder}/${pop}/MERGED/CLEANED/RECODED/chr${chr}.gen_info

		zcat ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen.gz | awk '{
		if (length($4)!=length($5)){
		if(length($4)>length($5)){
		printf $1" "$2" "$3" R D ";for(i=6;i<=NF;i++) printf "%s ", $i; printf "\n"
		} else {
		printf $1" "$2" "$3" R I ";for(i=6;i<=NF;i++) printf "%s ", $i; printf "\n"
		}
		} else {
		print $0
		}}' | gzip -c > ${basefolder}/${pop}/MERGED/CLEANED/RECODED/chr${chr}.gen.gz
		;;
	STEP4 )
	echo "Fourth step: Bimbam conversion"
		;;
	STEP5 )
	echo "Fifth step: Filevector conversion"
		;;
	CLEAN )
	echo "Clean temporary files"
		;;
	
esac




