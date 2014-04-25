#!/usr/local/bin/bash

#sed '26i\##FORMAT=<ID=PQ,Number=1,Type=Float,Description="IMPUTE2 phasing confidence score">' /nfs/team151_data03/VBSEQ/REL-2012-03-20/v1/esgi-vbseq.beagle.impute2.anno.20120427.csq.vcf > /lustre/scratch103/sanger/mc14/GENOTIPI/SEQUENCED/SEQUENCES/esgi-vbseq.beagle.impute2.anno.20120427.csq.vcf
#sed -i '27i\##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">' /lustre/scratch103/sanger/mc14/GENOTIPI/SEQUENCED/SEQUENCES/esgi-vbseq.beagle.impute2.anno.20120427.csq.vcf
sed -i '28i\##reference=file:///lustre/scratch105/projects/g1k/ref/main_project/human_g1k_v37.fasta' /lustre/scratch103/sanger/mc14/GENOTIPI/SEQUENCED/SEQUENCES/esgi-vbseq.beagle.impute2.anno.20120427.csq.vcf
