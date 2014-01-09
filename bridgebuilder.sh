#!/bin/bash

###############
#BAKER STAGE  #
###############
#Build baker single-ended way to avoid resorting the original BAM according to read names
if [ ! -f ${1}\_bridge.bam ]; then
		/software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa aln -q 15 -b /lustre/scratch113/projects/fvg_seq/F12HPCEUHK0358_HUMpngR/BRIDGED/makebridge_bridge.fa ${1} > ${1}_bridge.sai
		/software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa samse /lustre/scratch113/projects/fvg_seq/F12HPCEUHK0358_HUMpngR/BRIDGED/makebridge_bridge.fa ${1}_bridge.sai ${1} |samtools view -h -F 4 -Sb - > ${1}_bridge.bam
        rm ${1}\_bridge.sai
fi

###############
#Binnie STAGE #
###############
#2.1 Run binnie
echo "Run Binnie-ing"
~jr17/local/bin/binnie -i -s 50000 -m 10000 -v -a ${1} ${1}\_bridge.bam

#Run reheader
bam=${1}
remap=$(echo $bam"_remap.bam")
bridged=$(echo $bam"_bridged.bam")
unchanged=$(echo $bam"_unchanged.bam")

if [ -f ${remap} ]; then
        echo "Binnie-ing finished"
        #we need to sort remap.bam according to read name before realignment
        samtools sort -n $remap $remap\_readname\_sorted
        mv $remap\_readname\_sorted.bam $remap
        #2.2 remap reads in the *_remap.bam to target reference (hs37d5) pair-end
        # /software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa aln -q 15 -b1 /lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa/hs37d5.fa $remap > $remap.1.sai;
        /software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa aln -q 15 -b1 /lustre/scratch113/projects/crohns/2013jan04/bridge_builder_testing/bundle/hs37d5/hs37d5.fa $remap > $remap.1.sai;
        /software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa aln -q 15 -b2 /lustre/scratch113/projects/crohns/2013jan04/bridge_builder_testing/bundle/hs37d5/hs37d5.fa $remap > $remap.2.sai;
        # /software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa aln -q 15 -b2 /lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa/hs37d5.fa $remap > $remap.2.sai;
        /software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa sampe /lustre/scratch113/projects/crohns/2013jan04/bridge_builder_testing/bundle/hs37d5/hs37d5.fa $remap.1.sai $remap.2.sai $remap $remap | samtools view -h -Sb - > $remap\_hs37d5.bam
        # /software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa sampe /lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa/hs37d5.fa $remap.1.sai $remap.2.sai $remap $remap | samtools view -h -Sb - > $remap\_hs37d5.bam
        echo "remapping finished"
        rm $remap.[12].sai
        rm $remap
fi


if [[ -f ${bridged} && -f ${remap}\_hs37d5.bam ]];then
        #2.3 sort the bridged.bam and remaped remap.bam
        samtools sort -n $bridged $bridged\_readname
        samtools fixmate $bridged\_readname.bam $bridged\_fixmate.bam
        samtools sort $bridged\_fixmate.bam $bridged\_sorted
        samtools fixmate $remap\_hs37d5.bam $remap\_hs37d5.bam\_fixmate.bam
        #now sort by position
        samtools sort $remap\_hs37d5.bam\_fixmate.bam $remap\_remapsorted
        echo "sorting finished"
        rm $bridged
        rm $bridged\_readname.bam
        rm $bridged\_fixmate.bam
        #rm $remap\_hs37d5.bam
        rm $remap\_hs37d5.bam\_fixmate.bam
fi

###############
#BRUNEL STAGE #
###############
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



