#!/usr/bin/env bash
#
#run the merging of the panels
#we need a separate scripts singe we want to run it as a job array using the chunk file
start_pos=`sed -n "${LSB_JOBINDEX}p" $1 | awk '{print $1}'`
end_pos=`sed -n "${LSB_JOBINDEX}p" $1 | awk '{print $2}'`

pop1=$2
pop2=$3
chr=$4
gen_map=$5
hap_1=$6
hap_2=$7
leg_1=$8
leg_2=$9
buffer=${10}

chunk_n=$(printf "%03d" ${LSB_JOBINDEX})

out_ref=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.INGI_REF.${pop1}_${pop2}.${chunk_n}

if [[ -s ${out_ref}.hap.gz ]]; then
	echo "existing chunk! ${out_ref}, no gzip!!"
	# gzip ${out_ref}.hap
	# gzip ${out_ref}.legend

elif [[ -s ${out_ref}.hap ]]; then
	echo "existing chunk! ${out_ref}, gzip only!!"
	gzip ${out_ref}.hap
	gzip ${out_ref}.legend

else
	/nfs/team151/software/impute_v2.3.2_x86_64_static/impute2 -allow_large_regions -m ${gen_map} -h ${hap_1} ${hap_2} -l ${leg_1} ${leg_2} -k_hap 10000 10000 -merge_ref_panels -merge_ref_panels_output_ref ${out_ref} -merge_ref_panels_output_gen ${out_ref} -int ${start_pos} ${end_pos} -Ne 20000 -buffer ${buffer} -i ${out_ref}.info
	# /nfs/team151/software/impute_v2.3.2_x86_64_static/impute2 -allow_large_regions -m ${gen_map} -h ${hap_1} ${hap_2} -l ${leg_1} ${leg_2} -k_hap 2000 2000 -merge_ref_panels -merge_ref_panels_output_ref ${out_ref} -merge_ref_panels_output_gen ${out_ref} -int ${start_pos} ${end_pos} -Ne 20000 -buffer ${buffer} -i ${out_ref}.info
	gzip ${out_ref}.hap
	gzip ${out_ref}.legend
fi
# /nfs/team151/software/impute_v2.3.2_x86_64_static/impute2 -allow_large_regions -m ${gen_map} -h ${hap_1} ${hap_2} -l ${leg_1} ${leg_2} -k_hap 2000 2000 -merge_ref_panels -merge_ref_panels_output_gen ${out_ref} -int ${start_pos} ${end_pos} -Ne 20000 -buffer ${buffer} -i ${out_ref}.info
