#!/usr/local/bin/bash
#
#Script to launch fastqc on each fastaq file for each sample..

if [ $# -lt 2 ]
then
	echo -e "**********************\nWRONG ARGUMENT NUMBER!!!\n**********************"
        echo -e "USAGE:\n fastqc_launcher.sh <sample_file_path> <output file path> \n"
	exit 1
fi

sample_file_path=$1
outdir=$2
file_count=0

mkdir -p $outdir

for file in `ls $sample_file_path/*.fastq.gz`
do
	filecount=$[filecount+1]
	echo $filecount
	echo ${sample_file_path#*Sample_}
	qsub -N $filecount_${sample_file_path#*Sample_} -e $outdir/$filecount_fastqc.err -o $outdir/$filecount_fastqc.log -q workq -- /lustre1/tools/bin/fastqc $file --outdir=$outdir
done
