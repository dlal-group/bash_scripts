set -e

in_file=$1
out_file=$2
in_jobid=$3

ref="/lustre/scratch113/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa"

create_target_jobid=`echo "/software/jre1.7.0_25/bin/java -Xmx7g -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${ref} -I ${in_file} -o ${out_file}.intervals --known /lustre/scratch113/resources/variation/Homo_sapiens/grch37/gatk-bundle/2.5/Mills_and_1000G_gold_standard.indels.b37.vcf --known /lustre/scratch113/resources/variation/Homo_sapiens/grch37/gatk-bundle/2.5/1000G_phase1.indels.b37.vcf" | bsub -w "done($in_jobid)" -q basement -o log/create_target.%J.o -e log/create_target.%J.e -M 7900 -R 'select[mem>7900] rusage[mem=7900]' | sed -n 's/Job <\([0-9]*\)>.*/\1/p'`
realign_jobid=`echo "/software/jre1.7.0_25/bin/java -Xmx7g -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar -T IndelRealigner -R ${ref} -I ${in_file} -targetIntervals ${out_file}.intervals -known /lustre/scratch113/resources/variation/Homo_sapiens/grch37/gatk-bundle/2.5/Mills_and_1000G_gold_standard.indels.b37.vcf -known /lustre/scratch113/resources/variation/Homo_sapiens/grch37/gatk-bundle/2.5/1000G_phase1.indels.b37.vcf -o ${out_file} --consensusDeterminationModel USE_READS" | bsub -q basement -w "done($create_target_jobid)" -o log/realign.%J.o -e log/realign.%J.e -M 7900 -R 'select[mem>7900] rusage[mem=7900]' | sed -n 's/Job <\([0-9]*\)>.*/\1/p'`
echo $realign_jobid
