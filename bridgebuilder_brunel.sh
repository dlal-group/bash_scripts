#!/bin/bash

###############
#BRUNEL STAGE #
###############
#ARGS:
#
#$1=bam file
#$2=remap bam file
bam=${1}
bam_name=$(basename ${1})
remap=$(echo ${2}/$bam_name"_remap.bam")
bridged=$(echo ${2}/$bam_name"_bridged.bam")
unchanged=$(echo ${2}/$bam_name"_unchanged.bam")

echo $bam
echo $bam_name
echo ${2}
echo $remap
echo $bridged
echo $unchanged

if [ -f $remap\_remapsorted.bam ]; then
        # 3.1 Preparing header
        samtools view -H $remap\_remapsorted.bam|grep "^@SQ\|^@HD" >$2/$bam_name\_header.txt
        samtools view -H $bam|grep -v "^@SQ\|^@HD" >>$2/$bam_name\_header.txt
        echo -e "@PG\tID:bridgebuilder\tPN:bridgebuilder\tVN:0.3.3\tDS:BridgeBuilder used to realign sequences to a new reference\tCL:/lustre/scratch109/crohns/bridgebuilder/runBridgeBuilder.sh" $bam >>$2/$bam_name\_header.txt
        echo -e $RG >> $2/$bam_name\_header.txt

        output=$(basename $bam)
        #3.2 Correct header for bridged file
        samtools reheader $2/$bam_name\_header.txt ${bridged}\_sorted.bam |samtools view  -h - |samtools view -Sb - > ${bridged}\_reheader.bam
        mv ${bridged}\_reheader.bam ${bridged}\_sorted.bam
        
        #3.3 Merge three sorted bams
        samtools merge -h $2/$bam_name\_header.txt -cp $2/$bam_name\_new.bam $unchanged $remap\_remapsorted.bam $bridged\_sorted.bam 
        rm $unchanged
        rm $bridged\_sorted.bam
        rm $remap\_remapsorted.bam

        #Add in RG info in the header
        samtools reheader $2/$bam_name\_header.txt  $2/$bam_name\_new.bam |samtools view -h -|samtools view -Sb - > $2/$bam_name\_reheader.bam
        mv $2/$bam_name\_reheader.bam $2/$bam_name\_new.bam
        samtools index $2/$bam_name\_new.bam

        
        rm $2/$bam_name\_header.txt
        echo "Finished Bridge Building. Have a nice day!"        
        
fi
