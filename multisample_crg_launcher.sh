#!/bin/bash

#script to submit multicall jobs in batches with dependency
if [ $# -lt 1 ]
then
	echo "ATTENTION!!Missing argument!!!"
	echo "Usage: "
	echo "multisample_crg_launcher.sh <bam list file> <output_folder>"
	exit 1
fi

#launch jobs by chromosome based on how many region on each of them
mkdir -p $2/LOGS

#define output folder
OUTF=$2

for chr in {1..22} X
#for chr in 3 9 13 16 20
#for chr in 20
do
	mkdir -p $2/${chr}/.jobs
	mkdir -p $2/${chr}/LOGS
	log=$2/${chr}/LOGS
	#reg_file="REGIONS/nimblegen_plus50_chr${chr}_r*.bed"
	reg_file="/nfs/users/xe/ggirotto/multisample/REGIONS/nimblegen_plus50_chr${chr}.bed"
#	reg_size=`ls REGIONS/nimblegen_plus50_chr${chr}_r*.bed| wc -l`
	#reg=0
#	while [ $reg -lt "$reg_size" ] 
#	do
#		reg=$[reg+1]
		echo ${reg_file}
		#echo "echo \"bash $1/multisample_crg_call.sh ${chr} ${reg_file}\" | qsub -N \"chr${chr}_multicall\" -o \"${log}/chr${chr}_multicall.o\" -e \"${log}/chr${chr}_multicall.e\" -l h_rt=200:00:00 -l h_vmem=30Gb -cwd -q long" > $1/${chr}/.jobs/${chr}_r${reg}_job.sh
		#echo "echo \"bash /nfs/users/xe/ggirotto/multisample/scripts/multisample_crg_call.sh ${chr} ${reg_file}\" | qsub -N \"chr${chr}_multicall\" -o \"${log}/chr${chr}_multicall.o\" -e \"${log}/chr${chr}_multicall.e\" -l h_rt=200:00:00 -l virtual_free=20Gb -cwd -q long -pe smp 8" > $1/${chr}/.jobs/${chr}_job.sh
		echo "echo \"bash multisample_crg_call.sh ${chr} ${reg_file} $1 $2\" | qsub -N \"chr${chr}_multicall\" -o \"${log}/chr${chr}_multicall.o\" -e \"${log}/chr${chr}_multicall.e\" -l h_rt=200:00:00 -l virtual_free=40Gb -cwd -q xe-el6 -pe smp 8" > $2/${chr}/.jobs/${chr}_job.sh
		if [ ! -f $OUTF/${chr}.multisampleinitial.allregions.snps.done ]
		then
			echo ${OUTF}
			bash $2/${chr}/.jobs/${chr}_job.sh
		else
			echo "Chr ${chr} already successfully processed!"
		fi

done

#now we need to launch the concat jobs for each chr
#qsub -N "test_multicall_concat_chr${chr}" -o "${log}/test_multicall_concat_chr${chr}.$JOB_ID.o" -e "${log}/test_multicall_concat_chr${chr}.$JOB_ID.e"  -hold_jid "chr${chr}_multicall" -l h_rt=80:00:00 -l virtual_free=16Gb -cwd -q long /nfs/users/xe/ggirotto/multisample/scripts/multisample_crg_chr_concat.sh ${chr} $1
#	qsub -N "test_multicall_concat_chr${chr}" -o "${log}/test_multicall_concat_chr${chr}.o" -e "${log}/test_multicall_concat_chr${chr}.e" -l h_rt=80:00:00 -cwd -q long ./multisample_crg_chr_concat.sh ${chr}
#done

#last step: concat all chr together and apply VQSR
#qsub -N "test_multicall_VQSR" -o "LOGS/test_multicall_VQSR.$JOB_ID.o" -e "LOGS/test_multicall_VQSR.$JOB_ID.e" -hold_jid "test_multicall_concat_chr*" -l h_rt=80:00:00 -l virtual_free=16Gb -cwd -q long /nfs/users/xe/ggirotto/multisample/scripts/multisample_crg_vqsr.sh $1

#check if we have all chr called:
chr_num=`ls $OUTF/*.multisampleinitial.allregions.snps.done | wc -l`
   
if [ ${chr_num} -eq "23" ]
then
	#if we already have all the chr files we need no dependecies for vqsr job
        touch $OUTF/1.call.done
	echo "Launch VQSR step"
	qsub -N "multicall_VQSR" -o "$2/LOGS/multicall_VQSR.o" -e "$2/LOGS/multicall_VQSR.e" -l h_rt=80:00:00 -l virtual_free=16Gb -cwd -q long /nfs/users/xe/ggirotto/max/scripts/bash_scripts/multisample_crg_vqsr.sh $2
else
	echo "Launch VQSR step"
	qsub -N "multicall_VQSR" -o "$2/LOGS/multicall_VQSR.o" -e "$2/LOGS/multicall_VQSR.e" -hold_jid "chr*_multicall" -l h_rt=80:00:00 -l virtual_free=16Gb -cwd -q long /nfs/users/xe/ggirotto/max/scripts/bash_scripts/multisample_crg_vqsr.sh $2
fi


