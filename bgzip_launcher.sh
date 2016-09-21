#!/usr/local/bin/bash

#bgzip -c /nfs/team151_data03/VBSEQ/REL-2012-03-20/v1/esgi-vbseq.beagle.impute2.anno.20120427.csq.vcf > /lustre/scratch103/sanger/mc14/GENOTIPI/SEQUENCED/SEQUENCES/esgi-vbseq.beagle.impute2.anno.20120427.csq.vcf.gz

#bsub -J "bgzip_job" -o "%J.log"  -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
#-q basement \
bgzip -c /nfs/team151_data03/VBSEQ/REL-2012-03-20/v1/esgi-vbseq.beagle.impute2.anno.20120427.csq.vcf > /lustre/scratch103/sanger/mc14/GENOTIPI/SEQUENCED/SEQUENCES/esgi-vbseq.beagle.impute2.anno.20120427.csq.vcf.gz
#bgzip -d /lustre/scratch103/sanger/mc14/GENOTIPI/SEQUENCED/SEQUENCES/esgi-vbseq.beagle.impute2.anno.20120427.csq.vcf.gz 

