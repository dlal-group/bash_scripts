#!/bin/bash

#launcher for BridgeBuilder from Yang's script

#$1= log dir
#$2= bam file path

mkdir -p ${1}

#bams=$(ls -1 ${2}/*/Alignment_result/*.recal.bam)
bams=${2}
for file in $bams
do
echo "Working on: $file"
###############
#BAKER STAGE  #
###############
#Build baker single-ended way to avoid resorting the original BAM according to read names
#bsub -J "baker_bridge_builder_s$(basename ${file})" -o "${1}/%J_bridge_builder_baker_s$(basename ${file}).log" -e "${1}/%J_bridge_builder_baker_s$(basename ${file}).err" -M7000000 -R"select[mem>7000] rusage[mem=7000]" -q basement -- ~/Work/bash_scripts/bridgebuilder_baker.sh ${file}

###############
#Binnie STAGE #
###############
#2.1 Run binnie
#bsub -J "binnie_bridge_builder_s$(basename ${file})" -w "ended(baker_bridge_builder_s$(basename ${file}))" -o "${1}/%J_bridge_builder_binnie_s$(basename ${file}).log" -e "${1}/%J_bridge_builder_binnie_s$(basename ${file}).err" -M7000000 -R"select[mem>7000] rusage[mem=7000]" -q basement -- ~/Work/bash_scripts/bridgebuilder_binnie.sh ${file}


###############
#BRUNEL STAGE #
###############
#bsub -J "brunel_bridge_builder_s$(basename ${file})" -w "ended(binnie_bridge_builder_s$(basename ${file}))" -o "${1}/%J_bridge_builder_brunel_s$(basename ${file}).log" -e "${1}/%J_bridge_builder_brunel_s$(basename ${file}).err" -M7000000 -R"select[mem>7000] rusage[mem=7000]" -q basement -- ~/Work/bash_scripts/bridgebuilder_brunel.sh ${file}
bsub -J "brunel_bridge_builder_s$(basename ${file})" -o "${1}/%J_bridge_builder_brunel_s$(basename ${file}).log" -e "${1}/%J_bridge_builder_brunel_s$(basename ${file}).err" -M7000000 -R"select[mem>7000] rusage[mem=7000]" -q basement -- ~/Work/bash_scripts/bridgebuilder_brunel.sh ${file}

done
