#!/bin/bash

###############
#BRUNEL STAGE #
###############
#ARGS:
#
#$1=bam file
#$2=remap bam file
bam=${1}
remap=$(echo $bam"_remap.bam")
bridged=$(echo $bam"_bridged.bam")
unchanged=$(echo $bam"_unchanged.bam")

if [ -f $remap\_remapsorted.bam ]; then
        #3.1 Preparing header
        samtools view -H $remap\_remapsorted.bam|grep "^@SQ\|^@HD" >$bam\_header.txt
        samtools view -H $bam|grep -v "^@SQ\|^@HD" >>$bam\_header.txt
        echo -e "@PG\tID:bridgebuilder\tPN:bridgebuilder\tVN:0.3.3\tDS:BridgeBuilder used to realign sequences to a new reference\tCL:/lustre/scratch109/crohns/bridgebuilder/runBridgeBuilder.sh" $bam >>$bam\_header.txt

        output=$(basename $bam)
        #3.2 Merge three sorted bams 
        /lustre/scratch113/projects/crohns/2013jan04/bridge_builder_testing/bridgebuilder/brunel/src/brunel $bam\_header.txt $unchanged $remap\_remapsorted.bam $bridged\_sorted.bam $bam\_new.bam
        #reads=$(samtools view /lustre/scratch113/projects/crohns/uk10k-data/${output}|wc -l)
        #echo $output $reads>>/lustre/scratch113/projects/crohns/uk10k-bridgebuilder/new.stats
        # rm $unchanged
        # rm $bridged\_sorted.bam
        # rm $remap\_remapsorted.bam
        # rm $bam\_header.txt
        echo "Finished Bridge Building. Have a nice day!"        
fi
