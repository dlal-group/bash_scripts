#!/usr/local/bin/bash
#indexer script

#ARGS:
#$1=bam filepath
#$2=out dir

out_dir=$2
mkdir -p $out_dir/FASTQC
# files=(`ls $1/*.bam`)
# index=$[LSB_JOBINDEX - 1]
file=$1
#samtools stats
#(echo "ref_seq_name	seq_length	#mapped_reads	#unmapped_reads";samtools idxstats ${files[$index]}) > ${files[$index]}.stats

#extract some info from the bam
#(echo "MAPPING_QUAL   INSERT_SIZE";samtools view ${files[$index]} | cut -f 5,9) > ${files[$index]}.extract.stats

#fastqc stats
# fastqc ${files[$index]} --outdir=$out_dir/FASTQC
fastqc ${file} --outdir=$out_dir/FASTQC
#length=${files[@]}

#for ((index=0; index <= length; index++))
#do
#stats from bamcheck
# out_name=`basename ${files[$index]}`
out_name=`basename ${file}`
mkdir -p $out_dir/BAMCHECK_STATS/
echo "Processing $out_name"
# /software/vertres/codebase/scripts/bamcheck -c 1,50,1 -d ${files[$index]} > $out_dir/BAMCHECK_STATS/$out_name.bamchek.stats
/software/vertres/codebase/scripts/bamcheck -c 1,50,1 -d ${file} > $out_dir/BAMCHECK_STATS/$out_name.bamchek.stats
#done
