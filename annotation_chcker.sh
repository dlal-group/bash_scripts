#!/usr/local/bin/bash

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
	echo "CHR ${i}"
	#tabix /lustre/scratch107/projects/esgi-vbseq/REL-2012-06-07/v2/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz chr ${i} | fgrep -v "INDEL" | cut -f 1-3 | wc -l
	#tabix /lustre/scratch107/projects/esgi-vbseq/REL-2012-06-07/v1/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz chr ${i} | fgrep -v "INDEL" | cut -f 1-3 | wc -l
	#tabix /lustre/scratch107/projects/esgi-vbseq/REL-2012-06-07/v1/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz chr ${i} | fgrep -v "INDEL" | cut -f 1-3 | grep "\.$" | wc -l
	#tabix /lustre/scratch107/projects/esgi-vbseq/REL-2012-06-07/v1/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz chr ${i} | fgrep -v "INDEL" | cut -f 1-3 | grep "\.$" | head

	#tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.vcf.gz chr ${i} | fgrep -v "INDEL" | cut -f 1-3 | grep "\.$" | wc -l
	tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.vcf.gz chr ${i} | fgrep -v "INDEL" | grep "	\.	" | wc -l
	#tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.vcf.gz chr ${i} | fgrep -v "INDEL" | wc -l
	#tabix /lustre/scratch107/projects/esgi-vbseq/REL-2012-06-07/v2/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz chr ${i} | fgrep -v "INDEL" | cut -f 1-3 | grep "\.$" | head
	#grep "^${i}	" VBI.23chr.snps.vcf | cut -f 1-3 | grep "\.$" | wc -l
	#grep "^${i}	" VBI.23chr.snps.reann.vcf | cut -f 1-3 | grep "\.$" | wc -l
done

