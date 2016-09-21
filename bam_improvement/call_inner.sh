set -e

in_file=$1
out_file=$2

ref="/lustre/scratch113/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa"
dbsnp="/lustre/scratch113/resources/variation/Homo_sapiens/grch37/dbsnp_141.vcf.gz"
chrom=" -L 6"

/software/jre1.7.0_25/bin/java -Xmx28g -Djava.io.tmpdir=./tmp -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${ref} -nct 4 --dbsnp $dbsnp -ERC GVCF${chrom} -o ${out_file} -I ${in_file}  -variant_index_type LINEAR -variant_index_parameter 128000
