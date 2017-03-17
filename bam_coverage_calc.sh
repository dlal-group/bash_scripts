#!/usr/local/bin/bash

files=(`ls $1/*.bam`)
index=$[LSB_JOBINDEX - 1]

bam_file_path=$1
out_path=$2

mkdir -p ${out_path}

out_name=`basename ${files[$index]}`

#extract per base coverage from bam file, than calculate mean coverage (use R for a complete summary for each sample?)

samtools depth ${files[$index]} | gzip -c > ${out_path}/${out_name}.coverage.gz

#call R to create a complete summary
contig=`zcat ${out_path}/${out_name}.coverage.gz | cut -f 1 | sort | uniq`

for i in $contig
do	
	mkdir -p ${out_path}/CONTIG_${i}
	zcat ${out_path}/${out_name}.coverage.gz | grep "^${i}	" > ${out_path}/CONTIG_${i}/${out_name}_chr${i}.coverage
	/software/bin/R CMD BATCH '--args '${out_path}/CONTIG_${i}/${out_name}_chr${i}.coverage'' /nfs/users/nfs_m/mc14/Work/r_scripts/summarize_coverage.R ${out_path}/${out_name}_${i}.Rout
	rm ${out_path}/CHR${i}/${out_name}_chr${i}.coverage
done

