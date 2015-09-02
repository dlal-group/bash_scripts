#!/bin/bash

###############
#Binnie STAGE #
###############
#Args:
#$1:bam file
#$2:out dir
#
#2.1 Run binnie
bam_name=$(basename ${1})
bam=${1}
remap=$(echo ${2}/${bam_name}"_remap.bam")
bridged=$(echo ${2}/${bam_name}"_bridged.bam")
unchanged=$(echo ${2}/${bam_name}"_unchanged.bam")
module unload hgi/bwa/0.5.9
module load hgi/bwa/latest
# module add hgi/bwa/0.5.10

echo "Run Binnie-ing"
~jr17/local/bin/binnie -i -s 50000 -m 10000 -v -a -u ${unchanged} -b ${bridged} -r ${remap} ${1} ${2}/${bam_name}_bridge.bam

#Run reheader

echo $bam
echo $bam_name
echo ${2}
echo $remap
echo $bridged
echo $unchanged

if [ -f ${remap} ]; then
        echo "Binnie-ing finished"
        REF=/lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa
        #we need to sort remap.bam according to read name before realignment
        # /software/hgi/pkglocal/samtools-1.2/bin/samtools sort -n $remap $remap\_readname\_sorted
        # # /lustre/scratch113/projects/crohns/software/samtools-fixmate/samtools fixmate $remap\_readname\_sorted.bam $remap\_fixmate.bam
        # #582096.dedup.realn.recal.bam_new.bam_remap.bam_fixmate.bam
        # /software/hgi/pkglocal/samtools-1.2/bin/samtools fixmate $remap\_readname\_sorted.bam $remap\_fixmate.bam 
        # mv $remap\_fixmate.bam $remap
        # #2.2 remap reads in the *_remap.bam to target reference (hs37d5) pair-end
        # # /lustre/scratch113/projects/crohns/software/bwa.0.5.10_fixes/bwa aln -q 15 -b1 /lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa/hs37d5.fa $remap > $remap.1.sai;
        # bwa aln -q 15 -b1 $REF $remap > $remap.1.sai;
        # # /lustre/scratch113/projects/crohns/software/bwa.0.5.10_fixes/bwa aln -q 15 -b2 /lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa/hs37d5.fa $remap > $remap.2.sai;
        # bwa aln -q 15 -b2 $REF $remap > $remap.2.sai;
        # /lustre/scratch113/projects/crohns/software/bwa.0.5.10_fixes/bwa sampe /lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa/hs37d5.fa $remap.1.sai $remap.2.sai $remap $remap | samtools view -h -Sb - > $remap\_hs37d5.bam
        #####################################################################################################
        #2/09/2015 -> finally sort out the problem with duplicate removal! so this is a temporary fix!
        # bwa sampe $REF $remap.1.sai $remap.2.sai $remap $remap | /software/hgi/pkglocal/samtools-1.2/bin/samtools view -h -Sb - > $remap\_hs37d5.bam
        #map single end reads

        bwa samse $REF $remap.1.sai $remap | /software/hgi/pkglocal/samtools-1.2/bin/samtools view -h -Sb - > $remap\_hs37d5.1.bam
        bwa samse $REF $remap.2.sai $remap | /software/hgi/pkglocal/samtools-1.2/bin/samtools view -h -Sb - > $remap\_hs37d5.2.bam
        
        #sort single end reads
        /software/hgi/pkglocal/samtools-1.2/bin/samtools sort $remap\_hs37d5.1.bam $remap\_hs37d5.1.bam_remapsorted
        /software/hgi/pkglocal/samtools-1.2/bin/samtools sort $remap\_hs37d5.2.bam $remap\_hs37d5.2.bam_remapsorted

        #merge them back
        /software/hgi/pkglocal/samtools-1.2/bin/samtools merge -f -cp $remap\_hs37d5.bam $remap\_hs37d5.1.bam_remapsorted.bam $remap\_hs37d5.2.bam_remapsorted.bam

        #fix header
        /software/hgi/pkglocal/samtools-1.2/bin/samtools view -H $remap\_hs37d5.bam |grep -v "^@PG\|^@HD" > $remap\_hs37d5.bam_header.txt

        #rehead file
        /software/hgi/pkglocal/samtools-1.2/bin/samtools reheader $remap\_hs37d5.bam_header.txt $remap\_hs37d5.bam |/software/hgi/pkglocal/samtools-1.2/bin/samtools view  -h - |/software/hgi/pkglocal/samtools-1.2/bin/samtools view -Sb - > $remap\_hs37d5.bam_reheader.bam
        mv $remap\_hs37d5.bam_reheader.bam $remap\_hs37d5.bam
        echo "remapping finished"

        rm $remap\_hs37d5.[12].bam_remapsorted
        # rm $remap.[12].sai
        # rm $remap
fi



if [[ -f ${bridged} && -f ${remap}\_hs37d5.bam ]];then
        #2.3 sort the bridged.bam and remaped remap.bam
        #uncomment the following three lines, to have the complete procedure!
        # /software/hgi/pkglocal/samtools-1.2/bin/samtools sort -n $bridged $bridged\_readname
        # /software/hgi/pkglocal/samtools-1.2/bin/samtools fixmate $bridged\_readname.bam $bridged\_fixmate.bam
        # /software/hgi/pkglocal/samtools-1.2/bin/samtools sort $bridged\_fixmate.bam $bridged\_sorted
        /software/hgi/pkglocal/samtools-1.2/bin/samtools fixmate $remap\_hs37d5.bam $remap\_hs37d5.bam\_fixmate.bam
        #now sort by position
        /software/hgi/pkglocal/samtools-1.2/bin/samtools sort $remap\_hs37d5.bam\_fixmate.bam $remap\_remapsorted
        echo "sorting finished"
        # rm $bridged
        # rm $bridged\_readname.bam
        # rm $bridged\_fixmate.bam
        #rm $remap\_hs37d5.bam
        # rm $remap\_hs37d5.bam\_fixmate.bam
fi
