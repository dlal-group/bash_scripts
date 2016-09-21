set -e

module add hgi/samtools/1.2
module add hgi/bcftools/1.2

outdir="/nfs/users/nfs_m/mc14/Work/bash_scripts/bam_improvement/data"

for file in `cat manifest`; do
base=`basename $file .bam`
#mark dups
markdup_id=`./markdup.sh ${file} ${outdir}/${base}_mkdup.bam`
#realign
realign_id=`./realign.sh ${outdir}/${base}_mkdup.bam ${outdir}/${base}_realign.bam $downsample_id`
#+bqsr
bqsr_realign_id=`./bqsr.sh ${outdir}/${base}_realign.bam ${outdir}/${base}_realign_bqsr.bam $realign_id`
#call that
./call.sh ${outdir}/${base}_realign_bqsr.bam ${outdir}/${base}_realign_bqsr.vcf.gz $bqsr_realign_id
done
