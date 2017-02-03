#!/usr/bin/env bash
#
# set -e

pop=$1
chr=$2
basefolder=$3
# basefolder=/netapp02/data/imputation/INGI_TGP3/impute
mode=$4
basefolder2=$5
outfolder=$6

case $mode in
	STEP1 )
	echo "First step: MERGE and spurious text removal"
		# last_chunk=`ls ${basefolder}/${pop}/${chr}/*.gz| sort| tail -n1`
		# last_chunk_n=`echo ${last_chunk}| cut -f 2 -d "."`
		gen_list=${basefolder}/${pop}/${chr}/to_merge.list
		#first, merge all chunks for gen and info
		while read line
		# for chunk in `seq $last_chunk_n`
		do
		# chunkStr=`printf "%02d" $chunk`
		file_name=`basename ${line}`
		chunkStr=`echo ${file_name}|cut -f 2 -d "."`
		echo ${chr} ${chunkStr}

		cat ${basefolder}/${pop}/${chr}/chr${chr}.${chunkStr}.gen.gz >> ${basefolder}/${pop}/chr${chr}.gen_tmp1.gz
		cat ${basefolder}/${pop}/${chr}/chr${chr}.${chunkStr}.gen_info | sed 's,'"/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/CARL_FVG_VBI_TGP3_ALL/${chr}/${chr}.INGI_REF.CARL_FVG_VBI_TGP3_ALL.*.legend.gz:"',,g' >> ${basefolder}/${pop}/chr${chr}.gen_tmp1_info

		#create the COMPLETE merged imputed files
		done < <(cat ${gen_list})
		echo "Creating merged complete files...."
		mkdir -p ${basefolder}/${pop}/MERGED/ALL
		mkdir -p ${basefolder}/${pop}/MERGED/CLEANED
		mkdir -p ${basefolder}/${pop}/MERGED/CLEANED/RECODED

		zcat ${basefolder}/${pop}/chr${chr}.gen_tmp1.gz | sed 's,'"/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/CARL_FVG_VBI_TGP3_ALL/${chr}/${chr}.INGI_REF.CARL_FVG_VBI_TGP3_ALL.*.legend.gz:"',,g'| gzip -c > ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen.gz
		(echo "snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0";fgrep -v -w rs_id ${basefolder}/${pop}/chr${chr}.gen_tmp1_info) > ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen_info
		echo "ENDED!!"
		;;
	STEP2 )
		echo "Second step: Multiallelic removal"
		mkdir -p ${basefolder}/${pop}/MERGED/ALL
		mkdir -p ${basefolder}/${pop}/MERGED/CLEANED
		mkdir -p ${basefolder}/${pop}/MERGED/CLEANED/RECODED

		#extract a list of sites to keep:this sites will be uniq by position!
		/home/cocca/scripts/bash_scripts/clean_impute_files.py ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen ${chr} ${pop} ${basefolder}/${pop}/MERGED/CLEANED

		#sort the keeplist
		sort -g -k2,2 ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen.to_keep -o ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen.to_keep
		#create the rsId keeplist for qctool
		cut -f 1 -d " " ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen.to_keep > ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen_rsID.to_keep

		#now we have the list of stuff we want to keep, so we'll just grep what we need
		(echo "snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0";fgrep -f ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen.to_keep ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen_info) > ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen_info

		#here we have some problems..
		qctool -g ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen.gz -incl-rsids ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen_rsID.to_keep -omit-chromosome -og ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen.gz 
		#create sample file
		cut -f -6 -d " " ${basefolder}/${pop}/${chr}/chr${chr}.05.gen_samples > ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen_samples
		echo "ENDED!!!"
		;;
	STEP3 )
	echo "Third step: Long allele names recode"
		mkdir -p ${basefolder}/${pop}/MERGED/CLEANED
		mkdir -p ${basefolder}/${pop}/MERGED/CLEANED/RECODED

		#and last thing we need is to RECODE long alleles names
		awk -v chr=${chr} '{
		if (length($4)!=length($5)){
		if(length($4)>length($5)){
		printf $1" "chr":"$3"_R_D "$3" R D ";for(i=6;i<=NF;i++) printf "%s ", $i; printf "\n"
		} else {
		printf $1" "chr":"$3"_R_I "$3" R I ";for(i=6;i<=NF;i++) printf "%s ", $i; printf "\n"
		}
		} else {
		print $0
		}}' ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen_info > ${basefolder}/${pop}/MERGED/CLEANED/RECODED/chr${chr}.gen_info
		
		zcat ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen.gz | awk -v chr=${chr} '{
		if (length($4)!=length($5)){
		if(length($4)>length($5)){
		printf $1" "chr":"$3"_R_D "$3" R D ";for(i=6;i<=NF;i++) printf "%s ", $i; printf "\n"
		} else {
		printf $1" "chr":"$3"_R_I "$3" R I ";for(i=6;i<=NF;i++) printf "%s ", $i; printf "\n"
		}
		} else {
		print $0
		}}' | gzip -c > ${basefolder}/${pop}/MERGED/CLEANED/RECODED/chr${chr}.gen.gz
		echo "ENDED!!"
		;;
	STEP4 )
		echo "Fourth step: Filevector conversion"
		mkdir -p ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR
		mkdir -p ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/dose
		mkdir -p ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/prob
		#here we can convert the RECODED version of the CLEANED files (no duplicates by position)
		/home/cocca/scripts/bash_scripts/impute2mach_launcher.sh ${chr} ${basefolder}/${pop}/MERGED/CLEANED/RECODED/chr${chr}.gen.gz ${basefolder}/${pop}/MERGED/CLEANED/RECODED/chr${chr}.gen_info ${basefolder}/${pop}/MERGED/CLEANED/chr${chr}.gen_samples ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR
		#create MAP files (move it in filevector creation step!!)
		(echo "SNP Position A0 A1 Rsq";tail -n+2 ${basefolder}/${pop}/MERGED/CLEANED/RECODED/chr${chr}.gen_info | cut -f 2-5,7 -d " ") > ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/dose/${pop}_INGI_TGP3_chr${chr}.map
		cp ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/dose/${pop}_INGI_TGP3_chr${chr}.map ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/prob/${pop}_INGI_TGP3_chr${chr}.info

		echo "ENDED!!"
		;;
	STEP5 )
	echo "Fifth step: Bimbam conversion"
		#we conver the unfiltered and not recoded files, to bimbam format
		sel=$5
		mkdir -p ${basefolder}/${pop}/MERGED/ALL/BIMBAM
		/home/cocca/scripts/bash_scripts/impute2bimbam_inner.sh ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen.gz ${chr} ${basefolder}/${pop}/MERGED/ALL/BIMBAM $sel
		echo "ENDED!!"
		;;
	CLEAN )
	echo "Clean temporary files from all steps and format files!"
		rm -r ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/${chr}
		rm ${basefolder}/${pop}/chr${chr}.gen_tmp1.gz
		rm ${basefolder}/${pop}/chr${chr}.gen_tmp1_info
		mv ${basefolder}/${pop}/MERGED/ALL/BIMBAM/chr${chr}.gen.gz.bimbam.gz ${basefolder}/${pop}/MERGED/ALL/BIMBAM/chr${chr}.bimbam.gz
		mv ${basefolder}/${pop}/MERGED/ALL/BIMBAM/chr${chr}.gen.gz.pos ${basefolder}/${pop}/MERGED/ALL/BIMBAM/chr${chr}.bimbam.pos
		mkdir -p ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/dose
		mkdir -p ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/prob
		mv ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/${chr}.dose.* ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/dose/
		mv ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/${chr}.prob.* ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/prob/
		mv ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/dose/${chr}.dose.fvd ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/dose/${pop}_INGI_TGP3_chr${chr}.fvd
		mv ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/dose/${chr}.dose.fvi ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/dose/${pop}_INGI_TGP3_chr${chr}.fvi
		mv ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/prob/${chr}.prob.fvd ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/prob/${pop}_INGI_TGP3_chr${chr}.fvd
		mv ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/prob/${chr}.prob.fvi ${basefolder}/${pop}/MERGED/CLEANED/RECODED/FILEVECTOR/prob/${pop}_INGI_TGP3_chr${chr}.fvi
		
		echo "ENDED!!"
		;;
	MERGE )
		#create output folder
		mkdir -p ${outfolder}/${pop^^}/MERGED/ALL/
		#merge data
		# qctool -g ${basefolder}/${pop^^}/${pop}/MERGED/ALL/chr${chr}.gen.gz -s ${basefolder}/${pop^^}/${pop}/MERGED/ALL/${pop^^}.gen_samples -g ${basefolder2}/${pop^^}/MERGED/ALL/chr${chr}.gen.gz -s ${basefolder2}/${pop^^}/MERGED/ALL/${pop^^}.gen_samples -snp-match-fields position,alleles -match-alleles-to-cohort1 -og ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.joined.gen -os ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.joined.sample -omit-chromosome -sort
		# gzip ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.joined.gen
		# #extract stats
		# qctool -g ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.joined.gen.gz -s ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.joined.sample -sample-stats ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.joined.sample_stats -os ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.joined.new.sample
		# qctool -g ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.joined.gen.gz -snp-stats ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.joined.gen.snp_stats

		#generate stats with snptest also....
		# /netapp/nfs/softwares/snptest_v2.5.2/snptest_v2.5.2 -summary_stats_only -data ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.joined.gen.gz ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.sample -o ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.joined.gen.SNPTEST.snp_stats
		/netapp/nfs/softwares/snptest_v2.5.2/snptest_v2.5.2 -summary_stats_only -data ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.gen.gz ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.sample -o ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.joined.gen.SNPTEST.snp_stats
		#we need only rsid chromosome position alleleA alleleB info columns
		header="snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0"
		(echo "${header}"; egrep -v "^#" ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.joined.gen.SNPTEST.snp_stats | tail -n+2 | cut -f 1,2,4-6,9,14-16,18 -d " "  | awk '{print $0, ($8+($9*2))/($10*2)}'| awk '{print $1,$2,$3,$4,$5,$11,$6,0,0,0,0,0}') > ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.joined.gen.SNPTEST.info_a1_freq


		;;
	INFO )
		# Step to create info style files from qctools calculated stats
		# 1) define the file header
		# outfolder=$5
		basefolder3=$7
		info_head="snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0"
		mkdir -p ${outfolder}/${pop^^}/MERGED/ALL
		# to get the real number of gentyoped sites, we need to extract the genotyped sites from the relative info files, than we get a uniq list
		# and use that list to tag the sites in the new files because the only true genotyped sites will be the overlapping sites in the new merged result
		#<(awk '$9==2' ${basefolder}/${pop^^}/${pop}/MERGED/ALL/chr${chr}.gen_info | cut -f 2 -d " ") #> ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.genotyped_one.list
		#<(awk '$9==2' ${basefolder2}/${pop^^}/MERGED/ALL/chr${chr}.gen_info | cut -f 2 -d " ") # > ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.genotyped_second.list
		 # | sort | uniq > ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.genotyped_merge.list
		comm -12 <(awk '$9==2' ${basefolder}/${pop^^}/${pop}/MERGED/ALL/chr${chr}.gen_info | cut -f 2 -d " "|sort) <(awk '$9==2' ${basefolder2}/${pop^^}/MERGED/ALL/chr${chr}.gen_info | cut -f 2 -d " "| sort) > ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.genotyped_merge.list
		# 2) select the relevant columns from the stats file 
		# (echo ${info_head};tail -n+2 ${basefolder}/${pop^^}/MERGED/ALL/chr${chr}.gen.snp_stats |awk 'BEGIN{OFS=" "}{$15=sprintf("%.3f",$15)}{$19=sprintf("%.3f",$19)}{if($6==$7) print $1,$2,$4,$5,$6,$15,$19,0;else print $1,$2,$4,$5,$6,1-$15,$19,0}' |awk 'BEGIN{OFS=" "}{if($2~"rs") print $0,2,0,0,0;else print $0,0,0,0,0}') > ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.gen_info
		(tail -n+2 ${basefolder3}/${pop^^}/MERGED/ALL/chr${chr}.gen.snp_stats |awk 'BEGIN{OFS=" "}{$15=sprintf("%.3f",$15)}{$19=sprintf("%.3f",$19)}{if($6==$7) print $1,$2,$4,$5,$6,$15,$19,0;else print $1,$2,$4,$5,$6,1-$15,$19,0}') > ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.gen_info_tmp
		fgrep -w -f ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.genotyped_merge.list ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.gen_info_tmp | awk 'BEGIN{OFS=" "}{ print $0,2,0,0,0}' > ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.gen_info_tmp_genotyped
		fgrep -v -w -f ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.genotyped_merge.list ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.gen_info_tmp | awk 'BEGIN{OFS=" "}{ print $0,0,0,0,0}' > ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.gen_info_tmp_not_genotyped
		(echo ${info_head};cat ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.gen_info_tmp_genotyped ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.gen_info_tmp_not_genotyped | sort -g -k3,3) > ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.gen_info

		mkdir -p ${outfolder}/${pop^^}/MERGED/ALL/TMP
		mv ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.genotyped_merge.list ${outfolder}/${pop^^}/MERGED/ALL/TMP/
		mv ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.gen_info_tmp ${outfolder}/${pop^^}/MERGED/ALL/TMP/
		mv ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.gen_info_tmp_genotyped ${outfolder}/${pop^^}/MERGED/ALL/TMP/
		mv ${outfolder}/${pop^^}/MERGED/ALL/chr${chr}.gen_info_tmp_not_genotyped ${outfolder}/${pop^^}/MERGED/ALL/TMP/
		;;
esac
