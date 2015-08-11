#!/usr/local/bin/bash
#using a job array approach 

# mkdir -p LOGS;size=`wc -l /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_after_BEAGLE.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_after_BEAGLE.list
# mkdir -p LOGS;size=`wc -l /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_raw_call.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_raw_call.list
# mkdir -p LOGS;size=`wc -l /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_after_SHAPEIT.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_after_SHAPEIT.list
# mkdir -p LOGS;size=`wc -l /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_VQSR_filt_call.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_VQSR_filt_call.list
# mkdir -p LOGS;size=`wc -l /lustre/scratch114/teams/soranzo/users/mc14/BP_WGS_checks/BP_complete/BP_complete.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch114/teams/soranzo/users/mc14/BP_WGS_checks/BP_complete/BP_complete.list
# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/carl_seq/07202015/variant_calling/pooled/all_check.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch113/projects/carl_seq/07202015/variant_calling/pooled/all_check.list

vcfchk_path=$1
outname=`basename ${vcfchk_path}`
i=`echo ${outname%.*}`
echo ${vcfchk_path}
echo ${outname}
echo ${i}

(echo -e "CHR\tVARIANTS\tSNPS\tSingleton_SNPS\tMONO_SNPS\tINDELS\tSingleton_INDELS\tMONO_INDELS\tMULTI_SITES\tMULTI_SNPS\tMULTI_INDELS";
# for i in {1..22} X
# do
snps=`grep "^SN" ${vcfchk_path} | grep "number of SNPs"|cut -f 4`
indels=`grep "^SN" ${vcfchk_path} | grep "number of indels"|cut -f 4`
sis_snps=`grep "^SiS" ${vcfchk_path}  | cut -f 4`
sis_indels=`grep "^SiS" ${vcfchk_path} | cut -f 7`;
multi_sites=`grep "^SN" ${vcfchk_path} | grep "number of multiallelic sites"|cut -f 4`;
multi_snps=`grep "^SN" ${vcfchk_path} | grep "number of multiallelic SNP sites"|cut -f 4`;
multi_indel=$[multi_sites - multi_snps]
tot_v=$[snps + indels]
mono_snp=`grep "^AF" ${vcfchk_path}| awk '$3==0' |cut -f 4`
mono_indels=`grep "^AF" ${vcfchk_path}| awk '$3==0' |cut -f 7`
tot_mono=$[mono_snp + mono_indels]
echo -e "${i}\t${tot_v}\t${snps}\t${sis_snps}\t${mono_snp}\t${indels}\t${sis_indels}\t${mono_indels}\t${multi_sites}\t${multi_snps}\t${multi_indel}";
# done
) > ${outname}.tab
