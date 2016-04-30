#!/usr/local/bin/bash
#using a job array approach 

# mkdir -p LOGS;size=`wc -l /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_after_BEAGLE.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_after_BEAGLE.list
# mkdir -p LOGS;size=`wc -l /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_raw_call.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_raw_call.list
# mkdir -p LOGS;size=`wc -l /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_after_SHAPEIT.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_after_SHAPEIT.list
# mkdir -p LOGS;size=`wc -l /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_VQSR_filt_call.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /nfs/users/nfs_m/mc14/lustre114_home/BP_WGS_checks/BP_VQSR_filt_call.list
# mkdir -p LOGS;size=`wc -l /lustre/scratch114/teams/soranzo/users/mc14/BP_WGS_checks/BP_complete/BP_complete.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch114/teams/soranzo/users/mc14/BP_WGS_checks/BP_complete/BP_complete.list
# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/carl_seq/07202015/variant_calling/pooled/all_check.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch113/projects/carl_seq/07202015/variant_calling/pooled/all_check.list
# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/carl_seq/05082015/ALL/variant_calling/pooled/vcfchk_no_filter.lists |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch113/projects/carl_seq/05082015/ALL/variant_calling/pooled/vcfchk_no_filter.lists
# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/carl_seq/05082015/ALL/variant_calling/pooled/vcfchk_vqsr_filt.lists |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch113/projects/carl_seq/05082015/ALL/variant_calling/pooled/vcfchk_vqsr_filt.lists
# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/carl_seq/variant_refinement/18092015_BEAGLE/vcfstats.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch113/projects/carl_seq/variant_refinement/18092015_BEAGLE/vcfstats.list
# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/carl_seq/variant_refinement/18092015_SHAPEIT/vcfstats.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch113/projects/carl_seq/variant_refinement/18092015_SHAPEIT/vcfstats.list
# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/esgi-vbseq/08092015/variant_calling/pooled/vcfstats.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch113/projects/esgi-vbseq/08092015/variant_calling/pooled/vcfstats.list

# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/TRIMMED/files.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/TRIMMED/files.list
# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL/TRIMMED/files.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL/TRIMMED/files.list
# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/carl_seq/variant_refinement/12112015_FILTERED_REL/TRIMMED/files.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch113/projects/carl_seq/variant_refinement/12112015_FILTERED_REL/TRIMMED/files.list
# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/TGP3/TSI/trim_files.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/TGP3/TSI/trim_files.list
# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/TGP3/TSI/TRIMMED/files.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/TGP3/TSI/TRIMMED/files.list
# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/CARL_FVG_VBI_TSI/05042016_ANNOTATE/files.list |cut -f 1 -d " "`; bsub -J "call_check[1-${size}]" -o "LOGS/%J_call_check.%I.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/vcfchk_fast_count.sh /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/CARL_FVG_VBI_TSI/05042016_ANNOTATE/files.list

vcf_path=$1
outname=`basename ${vcf_path}`
out_path=`dirname ${vcf_path}`
vcfchk_path=${out_path}/${outname}.vcfchk

i=`echo ${outname%.*}`
echo ${vcfchk_path}
echo ${outname}
echo ${i}

#lets calculate the stats here!!
bcftools stats -F /lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa -s - -v ${vcf_path} > ${vcfchk_path}

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
