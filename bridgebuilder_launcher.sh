#!/bin/bash

#launcher for BridgeBuilder from Yang's script

# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/fvg_seq/04092015/BRIDGED/new_batch_misaligned_files.list|cut -f 1 -d " "`;bsub -J "bridge[1-${size}]" -o "LOGS/%J_bridge.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bridgebuilder_launcher.sh /lustre/scratch113/projects/fvg_seq/04092015/BRIDGED/new_batch_misaligned_files.list /lustre/scratch113/projects/fvg_seq/04092015/BRIDGED
# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/fvg_seq/04092015/BRIDGED/new_batch_misaligned_files_1.list|cut -f 1 -d " "`;bsub -J "bridge[1-${size}]" -o "LOGS/%J_bridge.%I.o" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bridgebuilder_launcher.sh /lustre/scratch113/projects/fvg_seq/04092015/BRIDGED/new_batch_misaligned_files_1.list /lustre/scratch113/projects/carl_seq/05262015/BRIDGED_FVG
# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/fvg_seq/04092015/BRIDGED/new_batch_misaligned_files.list|cut -f 1 -d " "`;bsub -J "bridge[1-${size}]" -o "LOGS/%J_bridge.%I.o" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bridgebuilder_launcher.sh /lustre/scratch113/projects/fvg_seq/04092015/BRIDGED/new_batch_misaligned_files.list /lustre/scratch113/projects/carl_seq/05262015/BRIDGED_FVG

#$1= bam file path
#$2= out dir

file=$1
# file=/nfs/users/nfs_m/mc14/fvg_seq/NEWBATCH/591582/Alignment_result/591582.dedup.realn.recal.bam
echo $file
file_name=$(basename ${file})
out_dir=$2/${file_name%.*}
mkdir -p $out_dir

module unload hgi/bwa/0.7.10
# for file in $bams
# do
echo "Working on: $file"
###############
#BAKER STAGE  #
###############
#Build baker single-ended way to avoid resorting the original BAM according to read names
bsub -J "baker_bridge_builder_s$(basename ${file})" -o "${out_dir}/%J_bridge_builder_baker_s$(basename ${file}).log" -e "$out_dir/%J_bridge_builder_baker_s$(basename ${file}).err" -M9000 -R"select[mem>9000] rusage[mem=9000]" -q long -- ~/Work/bash_scripts/bridgebuilder_baker.sh ${file} ${out_dir}

###############
#Binnie STAGE #
###############
#2.1 Run binnie
bsub -J "binnie_bridge_builder_s$(basename ${file})" -w "ended(baker_bridge_builder_s$(basename ${file}))" -o "${out_dir}/%J_bridge_builder_binnie_s$(basename ${file}).log" -e "${out_dir}/%J_bridge_builder_binnie_s$(basename ${file}).err" -M9000 -R"select[mem>9000] rusage[mem=9000]" -q basement -- ~/Work/bash_scripts/bridgebuilder_binnie.sh ${file} ${out_dir}
# bsub -J "binnie_bridge_builder_s$(basename ${file})" -o "${out_dir}/%J_bridge_builder_binnie_s$(basename ${file}).log" -e "${out_dir}/%J_bridge_builder_binnie_s$(basename ${file}).err" -M9000 -R"select[mem>9000] rusage[mem=9000]" -q basement -- ~/Work/bash_scripts/bridgebuilder_binnie.sh ${file} ${out_dir}


###############
#BRUNEL STAGE #
###############
bsub -J "brunel_bridge_builder_s$(basename ${file})" -w "ended(binnie_bridge_builder_s$(basename ${file}))" -o "${out_dir}/%J_bridge_builder_brunel_s$(basename ${file}).log" -e "${out_dir}/%J_bridge_builder_brunel_s$(basename ${file}).err" -M9000 -R"select[mem>9000] rusage[mem=9000]" -q basement -- ~/Work/bash_scripts/bridgebuilder_brunel.sh ${file} ${out_dir}
# bsub -J "brunel_bridge_builder_s$(basename ${file})" -o "${out_dir}/%J_bridge_builder_brunel_s$(basename ${file}).log" -e "${out_dir}/%J_bridge_builder_brunel_s$(basename ${file}).err" -M7000000 -R"select[mem>7000] rusage[mem=7000]" -q basement -- ~/Work/bash_scripts/bridgebuilder_brunel.sh ${file}

# done
