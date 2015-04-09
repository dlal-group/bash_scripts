#!/usr/local/bin/bash

#PIPELINE SCIRPT FOR ALIGNMENT CHECK
#select what kind of check you want to perform
#retrieve the MODE parameter to select the correct operation

#mkdir -p LOGS;size=`wc -l new_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_check[1-${size}]" -o "LOGS/%J_bam_check.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh new_batch_files.list /lustre/scratch113/projects/carl_seq/04072015/ref_check CHECKREF 38
#mkdir -p LOGS;size=`wc -l new_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_check[1-${size}]" -o "LOGS/%J_bam_check.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh new_batch_files.list /lustre/scratch113/projects/carl_seq/04072015/ref_check CHECKREF2 /lustre/scratch113/projects/carl_seq/release/20140325/sample_improved_bams_hgi_2/SC_CARLSEQ5554942.bam

#mkdir -p LOGS;size=`wc -l old_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_check[1-${size}]" -o "LOGS/%J_bam_check.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh old_batch_files.list /lustre/scratch113/projects/carl_seq/04072015/ref_check CHECKREF 38
#mkdir -p LOGS;size=`wc -l old_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_check[1-${size}]" -o "LOGS/%J_bam_check.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh old_batch_files.list /lustre/scratch113/projects/carl_seq/04072015/ref_check CHECKREF2 /lustre/scratch113/projects/carl_seq/NEWBATCH/603/result_alignment/603.rmdup.bam

#extract info per sample from bamcheck results
#Args
#$1 = files path
#$2 = out dir
out_name=`basename $1`
file=$1
out_dir=$2
MODE=$3


echo ${out_name}
echo ${file}
echo ${out_dir}
echo ${MODE}

# output folder relative to the chromosome, if specified
mkdir -p LOGS
mkdir -p ${out_dir}


case $MODE in
	CHECKREF)
	#call GATK with a depth check command to check if there is an error with the reference selected
	module add hgi/gatk-protected/latest
	ref=$4
	case $ref in
		37)
		REF=/lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta
		;;
		37d5)
		REF=/lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa
		;;
		38)
		# REF=/lustre/scratch114/resources/ref/Homo_sapiens/GRCh38_15/Homo_sapiens.GRCh38_15.fa #this should give an error
		REF=/lustre/scratch113/projects/carl_seq/04072015/ref_check/Homo_sapiens.GRCh38_15.fa #this should give an error
		;;
	esac
	
	#/software/hgi/pkglocal/gatk-protected-3.3/bin/gatk-queue

	/software/jre1.7.0_25/bin/java -Xmx1000m -Xms1000m -server -XX:+UseSerialGC -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar -T DepthOfCoverage \
	-I ${file} \
	-R ${REF} \
	--omitDepthOutputAtEachBase \
	--omitLocusTable \
	-o ${out_dir}/${out_name}.coverage
    ;;
    
	CHECKREF2)
	file2=$4

	echo -e "comparing file headers through diff\nfile1:${file}\nfile2:${file2}"
	samtools view -H ${file}| grep "^@SQ"| cut -f 1-3 | diff - <(samtools view -H ${file2}| grep "^@SQ"| cut -f 1-3) > ${out_dir}/${out_name}.diff
	# samtools view -H /lustre/scratch113/projects/carl_seq/release/20140325/sample_improved_bams_hgi_2/SC_CARLSEQ5554942.bam| grep "^@SQ"| cut -f 1-3 | diff - <(samtools view -H /lustre/scratch113/projects/carl_seq/NEWBATCH/603/result_alignment/603.rmdup.bam| grep "^@SQ"| cut -f 1-3)



	;;

	BAMSTATS )
	# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
	out_name=`basename $1`
	out_dir=$2

	mkdir -p $out_dir/BAMCHECK_STATS/
	mkdir -p $out_dir/BAMCHECK_STATS/PLOTS/QC_GRIND/

	echo "Processing $out_name"
	#create bam stats if they're not present already
	# /software/vertres/codebase/scripts/bamcheck -c 1,50,1 -d $1 > $out_dir/BAMCHECK_STATS/$out_name.bamchek.stats
	/software/vertres/codebase/scripts/bamcheck $1 > $out_dir/BAMCHECK_STATS/$out_name.bamchek.stats

	#grep ^SN ${files[$index]} | cut -f 2- > $1/$out_name.summary

	# ACGT content per cycle. Use `grep ^GCC | cut -f 2-` to extract this part. The columns are: cycle, and A,C,G,T counts [%]
	#(echo "cycle	A_count	C_count	G_count	T_count";grep ^GCC ${files[$index]} | cut -f 2-) > $1/$out_name.atgc

	# Indel distribution. Use `grep ^ID | cut -f 2-` to extract this part. The columns are: length, number of insertions, number of deletions
	#(echo "length	insertion_n	deletion_n";grep ^ID ${files[$index]} | cut -f 2-) > $1/$out_name.indel_dist

	# Indels per cycle. Use `grep ^IC | cut -f 2-` to extract this part. The columns are: cycle, number of insertions (fwd), .. (rev) , number of deletions (fwd), .. (rev)
	#(echo "cycle	insertion_fwd_n	insertion_rev_n	deletion_fwd_n	deletion_rev_n";grep ^IC ${files[$index]} | cut -f 2-) > $1/$out_name.indel_per_cycle

	# Coverage distribution. Use `grep ^COV | cut -f 2-` to extract this part
	#grep ^COV ${files[$index]} | cut -f 2- > $1/$out_name.cov_dist

	#summarize stats and plot
	#/software/bin/R CMD BATCH '--args '$1/$out_name'' /nfs/users/nfs_m/mc14/Work/r_scripts/summarize_bamcheck_stats.r

	#plot-bamcheck -r /lustre/scratch111/resources/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta.gc_stats -p $1/PLOTS/QC_GRIND/${out_name%%.*}/ $1/${out_name}
	#WE NEED TO USE THE NEW ALIGNEMNT REFERENCE
	echo "Plotting $out_dir/BAMCHECK_STATS/$out_name.bamchek.stats"
	# plot-bamcheck -r /lustre/scratch111/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa.gc_stats -p $out_dir/BAMCHECK_STATS/PLOTS/QC_GRIND/${out_name%%.*}/ $out_dir/BAMCHECK_STATS/$out_name.bamchek.stats
	plot-bamcheck -p $out_dir/BAMCHECK_STATS/PLOTS/QC_GRIND/${out_name%%.*} $out_dir/BAMCHECK_STATS/$out_name.bamchek.stats
    pop_path=$3
    # in_path=$4
  ;;
esac


