set -e

in_file=$1
out_file=$2
in_jobid=$3

ref="/lustre/scratch113/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa"

calc_covar_jobid=`echo "/software/jre1.7.0_25/bin/java -Xmx7g -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${ref} -I $in_file -o ${out_file}.table -knownSites /lustre/scratch113/resources/variation/Homo_sapiens/grch37/dbsnp_141.vcf.gz" | bsub -q basement -w "done($in_jobid)" -o log/calc_covar.%J.o -e log/calc_covar.%J.e -M 7900 -R 'select[mem>7900] rusage[mem=7900]' | sed -n 's/Job <\([0-9]*\)>.*/\1/p'`
bqsr_jobid=`echo "/software/jre1.7.0_25/bin/java -Xmx7g -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar -T PrintReads -R ${ref} -I ${in_file} -BQSR ${out_file}.table -o ${out_file}" | bsub -q basement -w "done($calc_covar_jobid)" -o log/bqsr.%J.o -e log/bqsr.%J.e -M 7900 -R 'select[mem>7900] rusage[mem=7900]' | sed -n 's/Job <\([0-9]*\)>.*/\1/p'`
echo $bqsr_jobid
