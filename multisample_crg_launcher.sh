#!/bin/bash

#script to submit multicall jobs in batches with dependency
if [ $# -lt 4 ]
then
	echo "ATTENTION!!Missing argument!!!"
	echo "Usage: "
	echo "multisample_crg_launcher.sh -l <bam list file> -o <output_folder> -t <VARIANT TYPE> -c <CALLER> [-m <Output mode (for GATK)>| -f <filter option>]"
	echo "Use -f flag if you want to launch VQSR filtering after calling"
	exit 1
fi

while getopts l:o:t:c:m:f opt; do
	case $opt in
		l )
		#bam file list
		BAMLIST=$OPTARG
		echo "${OPTARG}"		
			;;
		o )
		#define output folder
		echo "${OPTARG}"		
		OUTF=$OPTARG
		#launch jobs by chromosome based on how many region on each of them
		mkdir -p $OUTF/LOGS
			;;
		t )
		#variant type
		TYPE=$OPTARG
		echo "${OPTARG}"		
			;;
		c )
		#select the caller
		CALLER=$OPTARG
		echo "${OPTARG}"		
			;;
		m )
		#define the output mode
		OMODE=$OPTARG
		echo "${OPTARG}"		
			;;
		f )
		#set a flag for vqsr filtering after call
		VQSR='true'
		echo "${OPTARG}"		
			;;
		x )
		#set a flag for vqsr filtering SKIPPING the variant calling: assumes varian calling succesfully performed 
		# you have to provide all path arguments and variant types and so on,to perform the filtering step
		VQSR='true'
		CALLING='false'
		echo "${OPTARG}"		
		;;
		chr)
		CHROM=$OPTARG
		;;
	esac

done

#we're going to perform the calling step only if required
if [ -z ${CALLING+x} ]
then
	echo "No calling step required!"
else
	#let's start with autosomal data
	if [-z ${CHROM+x}]
	then
		for chr in {1..22}
		do
			mkdir -p ${OUTF}/${chr}/.jobs
			mkdir -p ${OUTF}/${chr}/LOGS
			log=${OUTF}/${chr}/LOGS
			#reg_file="REGIONS/nimblegen_plus50_chr${chr}_r*.bed"
			#reg_filelist=`ls /nfs/users/xe/ggirotto/multisample/REGIONS/splitted/nimblegen_plus50_chr${chr}_*.bed`
			#reg_filelist=`ls /nfs/users/xe/ggirotto/multisample/REGIONS_150/splitted/exome_region_for_coverage_enriched_150bp_chr${chr}.bed`
			reg_filelist=`ls /nfs/users/xe/ggirotto/multisample/REGIONS_150/splitted/exome_region_for_coverage_enriched_150bp_chr${chr}_*.bed`

			for reg_file in $reg_filelist
			do
				echo ${reg_file}
				filename=`basename ${reg_file}`
				#reg=`echo $filename|cut -f 4 -d "_"|cut -f 1 -d "."`
				reg=`echo ${filename%.*}|awk 'BEGIN{FS="_"};{print $(NF)}'`
				
				#echo "echo \"bash $1/multisample_crg_call.sh ${chr} ${reg_file}\" | qsub -N \"chr${chr}_multicall\" -o \"${log}/chr${chr}_multicall.o\" -e \"${log}/chr${chr}_multicall.e\" -l h_rt=200:00:00 -l h_vmem=30Gb -cwd -q long" > $1/${chr}/.jobs/${chr}_r${reg}_job.sh
				#echo "echo \"bash /nfs/users/xe/ggirotto/multisample/scripts/multisample_crg_call.sh ${chr} ${reg_file}\" | qsub -N \"chr${chr}_multicall\" -o \"${log}/chr${chr}_multicall.o\" -e \"${log}/chr${chr}_multicall.e\" -l h_rt=200:00:00 -l virtual_free=20Gb -cwd -q long -pe smp 8" > $1/${chr}/.jobs/${chr}_job.sh
				#echo "echo \"bash multisample_crg_call.sh ${chr} ${reg_file} $1 ${OUTF} ${TYPE} ${CALLER}\" | qsub -N \"chr${chr}_multicall\" -o \"${log}/chr${chr}_multicall.o\" -e \"${log}/chr${chr}_multicall.e\" -l h_rt=200:00:00 -l vf=40G -cwd -q xe-el6 -pe smp 8" > ${OUTF}/${chr}/.jobs/${chr}_job.${TYPE}.sh
				case ${CALLER} in
					GATK)
					echo "echo \"bash multisample_crg_call_${CALLER}.sh ${chr} ${reg_file} ${BAMLIST} ${OUTF} ${TYPE} ${OMODE}\" | qsub -N \"chr${chr}_${reg}_${CALLER}_multicall\" -o \"${log}/chr${chr}_${reg}_${CALLER}_multicall.o\" -e \"${log}/chr${chr}_${reg}_${CALLER}_multicall.e\" -l h_rt=200:00:00 -l vf=30G -cwd -q xe-el6 -pe smp 4" > ${OUTF}/${chr}/.jobs/${chr}_${reg}_${CALLER}_job.${TYPE}.sh
					;;
					
					SAMTOOLS)
					echo "echo \"bash multisample_crg_call_${CALLER}.sh ${chr} ${reg_file} ${BAMLIST} ${OUTF} ${TYPE} \" | qsub -N \"chr${chr}_${reg}_${CALLER}_multicall\" -o \"${log}/chr${chr}_${reg}_${CALLER}_multicall.o\" -e \"${log}/chr${chr}_${reg}_${CALLER}_multicall.e\" -l h_rt=200:00:00 -l vf=10G -cwd -q xe-el6" > ${OUTF}/${chr}/.jobs/${chr}_${reg}_${CALLER}_job.${TYPE}.sh
					;;
				esac

				# if [ ! -f $OUTF/${chr}.multisampleinitial.allregions.${TYPE}.done ]
				if [ ! -s $OUTF/${chr}/${chr}.${reg}.multisampleinitial.allregions.${TYPE}.done ]
				then
					echo ${OUTF}
					bash ${OUTF}/${chr}/.jobs/${chr}_${reg}_${CALLER}_job.${TYPE}.sh
				else
					echo "Chr${chr} region ${reg} , caller ${CALLER} already successfully processed for ${TYPE}s!"
				fi
			done
			#now we need to launch the concat jobs for each chr
			echo "Launch CHR${chr} concat step"
			qsub -N "multicall_concat_chr${chr}_${TYPE}" -o "$2/LOGS/multicall_concat_chr${chr}_${TYPE}.o" -e "$2/LOGS/multicall_concat_chr${chr}_${TYPE}.e" -hold_jid "chr${chr}_*_${CALLER}_multicall" -l h_rt=80:00:00 -l vf=16G -cwd -q xe-el6 /nfs/users/xe/ggirotto/max/scripts/bash_scripts/multisample_crg_reg_concat.sh $OUTF/${chr}/${chr} ${TYPE}
			#qsub -N "test_multicall_concat_chr${chr}" -o "${log}/test_multicall_concat_chr${chr}.$JOB_ID.o" -e "${log}/test_multicall_concat_chr${chr}.$JOB_ID.e"  -hold_jid "chr${chr}_multicall" -l h_rt=80:00:00 -l virtual_free=16Gb -cwd -q long /nfs/users/xe/ggirotto/multisample/scripts/multisample_crg_chr_concat.sh ${chr} $1
			#qsub -N "test_multicall_concat_chr${chr}" -o "${log}/test_multicall_concat_chr${chr}.o" -e "${log}/test_multicall_concat_chr${chr}.e" -l h_rt=80:00:00 -cwd -q long ./multisample_crg_chr_concat.sh ${chr}
			#done
		done
	else
		echo "Launch CHR${CHROM} calling"
		chr=${CHROM}
		mkdir -p ${OUTF}/${chr}/.jobs
		mkdir -p ${OUTF}/${chr}/LOGS
		log=${OUTF}/${chr}/LOGS
		reg_filelist=`ls /nfs/users/xe/ggirotto/multisample/REGIONS_150/splitted/exome_region_for_coverage_enriched_150bp_chr${chr}_*.bed`

		for reg_file in $reg_filelist
		do
			echo ${reg_file}
			filename=`basename ${reg_file}`
			#reg=`echo $filename|cut -f 4 -d "_"|cut -f 1 -d "."`
			reg=`echo ${filename%.*}|awk 'BEGIN{FS="_"};{print $(NF)}'`
			
			#echo "echo \"bash $1/multisample_crg_call.sh ${chr} ${reg_file}\" | qsub -N \"chr${chr}_multicall\" -o \"${log}/chr${chr}_multicall.o\" -e \"${log}/chr${chr}_multicall.e\" -l h_rt=200:00:00 -l h_vmem=30Gb -cwd -q long" > $1/${chr}/.jobs/${chr}_r${reg}_job.sh
			#echo "echo \"bash /nfs/users/xe/ggirotto/multisample/scripts/multisample_crg_call.sh ${chr} ${reg_file}\" | qsub -N \"chr${chr}_multicall\" -o \"${log}/chr${chr}_multicall.o\" -e \"${log}/chr${chr}_multicall.e\" -l h_rt=200:00:00 -l virtual_free=20Gb -cwd -q long -pe smp 8" > $1/${chr}/.jobs/${chr}_job.sh
			#echo "echo \"bash multisample_crg_call.sh ${chr} ${reg_file} $1 ${OUTF} ${TYPE} ${CALLER}\" | qsub -N \"chr${chr}_multicall\" -o \"${log}/chr${chr}_multicall.o\" -e \"${log}/chr${chr}_multicall.e\" -l h_rt=200:00:00 -l vf=40G -cwd -q xe-el6 -pe smp 8" > ${OUTF}/${chr}/.jobs/${chr}_job.${TYPE}.sh
			case ${CALLER} in
				GATK)
				echo "echo \"bash multisample_crg_call_${CALLER}.sh ${chr} ${reg_file} ${BAMLIST} ${OUTF} ${TYPE} ${OMODE}\" | qsub -N \"chr${chr}_${reg}_${CALLER}_multicall\" -o \"${log}/chr${chr}_${reg}_${CALLER}_multicall.o\" -e \"${log}/chr${chr}_${reg}_${CALLER}_multicall.e\" -l h_rt=200:00:00 -l vf=30G -cwd -q xe-el6 -pe smp 4" > ${OUTF}/${chr}/.jobs/${chr}_${reg}_${CALLER}_job.${TYPE}.sh
				;;
				
				SAMTOOLS)
				echo "echo \"bash multisample_crg_call_${CALLER}.sh ${chr} ${reg_file} ${BAMLIST} ${OUTF} ${TYPE} \" | qsub -N \"chr${chr}_${reg}_${CALLER}_multicall\" -o \"${log}/chr${chr}_${reg}_${CALLER}_multicall.o\" -e \"${log}/chr${chr}_${reg}_${CALLER}_multicall.e\" -l h_rt=200:00:00 -l vf=10G -cwd -q xe-el6" > ${OUTF}/${chr}/.jobs/${chr}_${reg}_${CALLER}_job.${TYPE}.sh
				;;
			esac

			if [ ! -s $OUTF/${chr}/${chr}.${reg}.multisampleinitial.allregions.${TYPE}.done ]
			then
				echo ${OUTF}
				bash ${OUTF}/${chr}/.jobs/${chr}_${reg}_${CALLER}_job.${TYPE}.sh
			else
				echo "Chr${chr} region ${reg} , caller ${CALLER} already successfully processed for ${TYPE}s!"
			fi
		done
		#now we need to launch the concat jobs for each chr
		echo "Launch CHR${chr} concat step"
		qsub -N "multicall_concat_chr${chr}_${TYPE}" -o "$2/LOGS/multicall_concat_chr${chr}_${TYPE}.o" -e "$2/LOGS/multicall_concat_chr${chr}_${TYPE}.e" -hold_jid "chr${chr}_*_${CALLER}_multicall" -l h_rt=80:00:00 -l vf=16G -cwd -q xe-el6 /nfs/users/xe/ggirotto/max/scripts/bash_scripts/multisample_crg_reg_concat.sh $OUTF/${chr}/${chr} ${TYPE}
		#qsub -N "test_multicall_concat_chr${chr}" -o "${log}/test_multicall_concat_chr${chr}.$JOB_ID.o" -e "${log}/test_multicall_concat_chr${chr}.$JOB_ID.e"  -hold_jid "chr${chr}_multicall" -l h_rt=80:00:00 -l virtual_free=16Gb -cwd -q long /nfs/users/xe/ggirotto/multisample/scripts/multisample_crg_chr_concat.sh ${chr} $1
		#qsub -N "test_multicall_concat_chr${chr}" -o "${log}/test_multicall_concat_chr${chr}.o" -e "${log}/test_multicall_concat_chr${chr}.e" -l h_rt=80:00:00 -cwd -q long ./multisample_crg_chr_concat.sh ${chr}
	fi
fi


#last step: concat all chr together and apply VQSR
#qsub -N "test_multicall_VQSR" -o "LOGS/test_multicall_VQSR.$JOB_ID.o" -e "LOGS/test_multicall_VQSR.$JOB_ID.e" -hold_jid "test_multicall_concat_chr*" -l h_rt=80:00:00 -l virtual_free=16Gb -cwd -q long /nfs/users/xe/ggirotto/multisample/scripts/multisample_crg_vqsr.sh $1


if [ -z ${VQSR+x} ]
then
	echo "No filtering step required!"
else
	#check if we have all chr called and merged:if not, we'll have some errors and use job dependency
	echo "Filtering step required!"
	chr_num=`ls ${OUTF}/*.multisampleinitial.allregions.${TYPE}.done | wc -l`
	echo "Done ${chr_num} chr!"
	
	if [ ${chr_num} -eq "23" ]
	then
		#if we already have all the chr files we need no dependecies for vqsr job
		touch $OUTF/1.call.${TYPE}.done
		echo "Launch VQSR step"
		qsub -N "multicall_VQSR" -o "$2/LOGS/multicall_${TYPE}_VQSR.o" -e "$2/LOGS/multicall_${TYPE}_VQSR.e" -l h_rt=80:00:00 -l vf=16G -cwd -q xe-el6 /nfs/users/xe/ggirotto/max/scripts/bash_scripts/multisample_crg_vqsr.sh $2 ${TYPE}
	else
		echo "Launch VQSR step"
		#qsub -N "multicall_VQSR" -o "$2/LOGS/multicall_${TYPE}_VQSR.o" -e "$2/LOGS/multicall_${TYPE}_VQSR.e" -hold_jid "chr*_multicall" -l h_rt=80:00:00 -l vf=16G -cwd -q xe-el6 /nfs/users/xe/ggirotto/max/scripts/bash_scripts/multisample_crg_vqsr.sh $2 ${TYPE}
		qsub -N "multicall_VQSR" -o "$2/LOGS/multicall_${TYPE}_VQSR.o" -e "$2/LOGS/multicall_${TYPE}_VQSR.e" -hold_jid "multicall_concat_chr*_${TYPE}" -l h_rt=80:00:00 -l vf=16G -cwd -q xe-el6 /nfs/users/xe/ggirotto/max/scripts/bash_scripts/multisample_crg_vqsr.sh $2 ${TYPE}
	fi
	
fi

