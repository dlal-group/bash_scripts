#!/usr/bin/env bash
#
#generate chunks for panel merging based on current legends
function find_min(){
pos_array=($@)
a_min=${pos_array[0]}

# Loop through all elements in the array
for i in "${pos_array[@]}"
do
    # Update min if applicable
    if [[ "$i" -lt "$a_min" ]]; then
        a_min="$i"
    fi
done
echo ${a_min}

}

function find_max(){
pos_array=($@)
a_max=${pos_array[0]}

# Loop through all elements in the array
for i in "${pos_array[@]}"
do
    # Update max if applicable
    if [[ "$i" -gt "$a_max" ]]; then
        a_max="$i"
    fi
done

echo ${a_max}
}

###############################################################
leg_1=$1
leg_2=$2
chunk_size=$3
pop1=$4
pop2=$5
chr=$6

s1=`zcat ${leg_1} | head -2 | tail -n1 | awk '{print $2}'`
s2=`zcat ${leg_2} | head -2 | tail -n1 | awk '{print $2}'`
starts=(${s1} ${s2})

e1=`zcat ${leg_1} | tail -n1 | awk '{print $2}'`
e2=`zcat ${leg_2} | tail -n1 | awk '{print $2}'`
ends=(${e1} ${e2})

start_chr=`find_min "${starts[@]}"`
end_chr=`find_max "${ends[@]}"`

start_pos=${start_chr}
end_pos=0
# than create the chunk file to use to submit the job array
while [ ${start_pos} -lt ${end_chr} ]
do
end_pos=$[start_pos + chunk_size]
if [ ${end_pos} -ge ${end_chr} ];then
end_pos=${end_chr}
fi
echo "${start_pos} ${end_pos}"
start_pos=$[end_pos + 1 ] 
done > /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.chunks.txt