#!/bin/bash

###############
#Binnie STAGE #
###############
#Args:
#$1:bam file
#
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
        /lustre/scratch113/projects/crohns/software/samtools-fixmate/samtools fixmate $remap\_readname\_sorted.bam $remap\_fixmate.bam
        mv $remap\_fixmate.bam $remap
        #2.2 remap reads in the *_remap.bam to target reference (hs37d5) pair-end
        /lustre/scratch113/projects/crohns/software/bwa.0.5.10_fixes/bwa aln -q 15 -b1 /lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa/hs37d5.fa $remap > $remap.1.sai;
        /lustre/scratch113/projects/crohns/software/bwa.0.5.10_fixes/bwa aln -q 15 -b2 /lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa/hs37d5.fa $remap > $remap.2.sai;
        /lustre/scratch113/projects/crohns/software/bwa.0.5.10_fixes/bwa sampe /lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa/hs37d5.fa $remap.1.sai $remap.2.sai $remap $remap | samtools view -h -Sb - > $remap\_hs37d5.bam
        echo "remapping finished"
        rm $remap.[12].sai
        rm $remap
fi


# if [ -f ${remap} ]; then
#         echo "Binnie-ing finished"
#         #we need to sort remap.bam according to read name before realignment
#         samtools sort -n $remap $remap\_readname\_sorted
#         mv $remap\_readname\_sorted.bam $remap
#         #2.2 remap reads in the *_remap.bam to target reference (hs37d5) pair-end
#         # /software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa aln -q 15 -b1 /lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa/hs37d5.fa $remap > $remap.1.sai;
#         /software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa aln -q 15 -b1 /lustre/scratch113/projects/crohns/2013jan04/bridge_builder_testing/bundle/hs37d5/hs37d5.fa $remap > $remap.1.sai;
#         /software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa aln -q 15 -b2 /lustre/scratch113/projects/crohns/2013jan04/bridge_builder_testing/bundle/hs37d5/hs37d5.fa $remap > $remap.2.sai;
#         # /software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa aln -q 15 -b2 /lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa/hs37d5.fa $remap > $remap.2.sai;
#         /software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa sampe /lustre/scratch113/projects/crohns/2013jan04/bridge_builder_testing/bundle/hs37d5/hs37d5.fa $remap.1.sai $remap.2.sai $remap $remap | samtools view -h -Sb - > $remap\_hs37d5.bam
#         # /software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa sampe /lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa/hs37d5.fa $remap.1.sai $remap.2.sai $remap $remap | samtools view -h -Sb - > $remap\_hs37d5.bam
#         echo "remapping finished"
#         rm $remap.[12].sai
#         rm $remap
# fi


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
