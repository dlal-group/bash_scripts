#!/usr/local/bin/bash

#PIPELINE SCIRPT FOR ALIGNMENT CHECK
#select what kind of check you want to perform
#retrieve the MODE parameter to select the correct operation
#add HGI modules
module add hgi/bcftools/latest
module add hgi/samtools/latest

# mkdir -p LOGS;size=`wc -l new_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_check[1-${size}]" -o "LOGS/%J_bam_check.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh new_batch_files.list /lustre/scratch113/projects/carl_seq/04072015/ref_check CHECKREF 38
# mkdir -p LOGS;size=`wc -l new_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_check[1-${size}]" -o "LOGS/%J_bam_check.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh new_batch_files.list /lustre/scratch113/projects/carl_seq/04072015/ref_check CHECKREF2 /lustre/scratch113/projects/carl_seq/release/20140325/sample_improved_bams_hgi_2/SC_CARLSEQ5554942.bam
# mkdir -p LOGS;size=`wc -l new_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_check[1-${size}]" -o "LOGS/%J_bam_check.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh new_batch_files.list /lustre/scratch113/projects/carl_seq/04072015/fofn FOFN carl_seq

# mkdir -p LOGS;size=`wc -l new_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_check[1-${size}]" -o "LOGS/%J_bam_check.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh new_batch_files.list /lustre/scratch113/projects/fvg_seq/04092015/ref_check CHECKREF2 /lustre/scratch113/projects/fvg_seq/release/20130827/sample_improved_bams_hgi/SC_FVGSEQ5316418.bam
# mkdir -p LOGS;size=`wc -l all_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_check[1-${size}]" -o "LOGS/%J_bam_check.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh all_batch_files.list /lustre/scratch113/projects/fvg_seq/04092015/bamcheck BAMSTATS
# mkdir -p LOGS;size=`wc -l all_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_check[1-${size}]" -o "LOGS/%J_bam_check.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh all_batch_files.list /lustre/scratch113/projects/carl_seq/04072015/bamcheck BAMSTATS


#mkdir -p LOGS;size=`wc -l old_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_check[1-${size}]" -o "LOGS/%J_bam_check.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh old_batch_files.list /lustre/scratch113/projects/carl_seq/04072015/ref_check CHECKREF 38
#mkdir -p LOGS;size=`wc -l old_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_check[1-${size}]" -o "LOGS/%J_bam_check.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh old_batch_files.list /lustre/scratch113/projects/carl_seq/04072015/ref_check CHECKREF2 /lustre/scratch113/projects/carl_seq/NEWBATCH/603/result_alignment/603.rmdup.bam
#mkdir -p LOGS;size=`wc -l old_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_check[1-${size}]" -o "LOGS/%J_bam_check.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh old_batch_files.list /lustre/scratch113/projects/fvg_seq/04092015/ref_check CHECKREF2 /lustre/scratch113/projects/fvg_seq/NEWBATCH/605889/result_alignment/605889.rmdup.bam


#extract info per sample from bamcheck results
#Args
#$1 = files path
#$2 = out dir
out_name=`basename $1`
file=$1
out_dir=$2
MODE=$3

echo "Parameters check before entering in the right mode:"
echo "Out name:${out_name}"
echo "File name:${file}"
echo "Out dir name:${out_dir}"
echo "MODE name:${MODE}"
echo "#######################################"
# output folder relative to the chromosome, if specified
# mkdir -p LOGS



case $MODE in
	IMPROVEDCHECK )
	# mkdir -p LOGS;size=`wc -l new_batch_files.list|cut -f 1 -d " "`;bsub -J "IMPROVEDCHECK[1-${size}]" -o "LOGS/%J_IMPROVEDCHECK.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh new_batch_files.list /lustre/scratch113/projects/carl_seq/05012015/improve_check IMPROVEDCHECK
	# mkdir -p LOGS;size=`wc -l new_batch_files.list|cut -f 1 -d " "`;bsub -J "IMPROVEDCHECK[1-${size}]" -o "LOGS/%J_IMPROVEDCHECK.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh new_batch_files.list /lustre/scratch113/projects/fvg_seq/05012015/improve_check IMPROVEDCHECK
	# mkdir -p LOGS;size=`wc -l new_batch_files.list|cut -f 1 -d " "`;bsub -J "IMPROVEDCHECK[1-${size}]" -o "LOGS/%J_IMPROVEDCHECK.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh new_batch_files.list /lustre/scratch113/projects/esgi-vbseq/05012015/improve_check IMPROVEDCHECK
	#extract PG tag from each file to check if they already underwent any improvement already
	mkdir -p ${out_dir}
	samtools view -H ${file} | grep "^@PG" > ${out_dir}/${out_name}.pgtag
	;;
	BAM2CRAMCHECK )
	# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/carl_seq/04072015/CRAMMED/CHECK/crammed_files.list|cut -f 1 -d " "`;bsub -J "cram_check[1-${size}]" -o "LOGS/%J_cram_check.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -d ~/Work/bash_scripts/bam_check_pipeline.sh /lustre/scratch113/projects/carl_seq/04072015/CRAMMED/CHECK/crammed_files.list BAM2CRAMCHECK /lustre/scratch113/projects/carl_seq/04072015/CRAMMED/CHECK
	
	# we should check if we've lost any data during the conversion from bam to cram
	# We need to compare: CIGAR strings, READs strings, QUALITY strings
	# To maximize parallelization we should submit a total of 6 jobs, extracting info from each file,
	# than calculating MD5sum for the column than compare the values and provide a final report
	file1=$1
	file2=$2
	out_dir=$4
	out_name=`basename ${file2}`
	mkdir -p ${out_dir}
	echo "Parameters check in ${MODE} :"
	echo "file1 : $file1"
	echo "file2 : $file2"
	echo "out_dir : $out_dir"
	echo "out_name : $out_name"
	echo "#######################################"
	# first check if we are comparing the right file

	#checks for bam file
	echo "samtools view ${file1}| cut -f 6 | md5sum | cut -f 1 -d ' ' > ${out_dir}/${out_name}.bam.CIGAR.md5" | bsub -J "check_${out_name}_CIGAR_bam" -o "LOGS/%J_check_${out_name}_CIGAR_bam.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal
	echo "samtools view ${file1}| cut -f 10 | md5sum | cut -f 1 -d ' ' > ${out_dir}/${out_name}.bam.READS.md5" | bsub -J "check_${out_name}_READS_bam" -o "LOGS/%J_check_${out_name}_READS_bam.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal
	echo "samtools view ${file1}| cut -f 11 | md5sum | cut -f 1 -d ' ' > ${out_dir}/${out_name}.bam.QUALITY.md5" | bsub -J "check_${out_name}_QUALITY_bam" -o "LOGS/%J_check_${out_name}_QUALITY_bam.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal

	#checks for cram file
	echo "samtools view ${file2}| cut -f 6 | md5sum | cut -f 1 -d ' ' > ${out_dir}/${out_name}.cram.CIGAR.md5" | bsub -J "check_${out_name}_CIGAR_cram" -o "LOGS/%J_check_${out_name}_CIGAR_cram.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal
	echo "samtools view ${file2}| cut -f 10 | md5sum | cut -f 1 -d ' ' > ${out_dir}/${out_name}.cram.READS.md5" | bsub -J "check_${out_name}_READS_cram" -o "LOGS/%J_check_${out_name}_READS_cram.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal
	echo "samtools view ${file2}| cut -f 11 | md5sum | cut -f 1 -d ' ' > ${out_dir}/${out_name}.cram.QUALITY.md5" | bsub -J "check_${out_name}_QUALITY_cram" -o "LOGS/%J_check_${out_name}_QUALITY_cram.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal
		
	echo "################ DIFF CHECK ###################"
	echo -e "diff_CIGAR=\`diff ${out_dir}/${out_name}.bam.CIGAR.md5 ${out_dir}/${out_name}.cram.CIGAR.md5\`\n
	diff_READS=\`diff ${out_dir}/${out_name}.bam.READS.md5 ${out_dir}/${out_name}.cram.READS.md5\`\n
	diff_QUALITY=\`diff ${out_dir}/${out_name}.bam.QUALITY.md5 ${out_dir}/${out_name}.cram.QUALITY.md5\`\n

	if [[ -z \"\${diff_CIGAR}\" ]]; then\n
		echo \"No differences after conversion between CIGAR data!! TOP!!!\" >> ${out_dir}/${out_name}.diff_report\n
	else\n
		echo \"D'oh!There are differences after conversion!! check your md5 CIGAR file!!!\" >> ${out_dir}/${out_name}.diff_report\n
	fi\n
	if [[ -z \"\${diff_READS}\" ]]; then\n
		echo \"No differences after conversion between READS data!! TOP!!!\" >> ${out_dir}/${out_name}.diff_report\n
	else\n
		echo \"D'oh!There are differences after conversion!! check your md5 READS file!!!\" >> ${out_dir}/${out_name}.diff_report\n
	fi	\n
	if [[ -z \"\${diff_QUALITY}\" ]]; then\n
		echo \"No differences after conversion between QUALITY data!! TOP!!!\" >> ${out_dir}/${out_name}.diff_report\n
	else\n
		echo \"D'oh!There are differences after conversion!! check your md5 QUALITY file!!!\" >> ${out_dir}/${out_name}.diff_report\n
	fi\n
	" | bsub -J "diff_check_${out_name}" -o "LOGS/%J_diff_check_${out_name}.o" -w "ended(check_${out_name}_*)" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal

	;;
	BAM2CRAM )
	# convert bam files to cram format to save space
	# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/carl_seq/04072015/CRAMMED/to_cram_file_1.list|cut -f 1 -d " "`;bsub -J "bam2cram[1-${size}]" -o "LOGS/%J_bam2cram.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh /lustre/scratch113/projects/carl_seq/04072015/CRAMMED/to_cram_file_1.list /lustre/scratch113/projects/carl_seq/04072015/CRAMMED/37d5 BAM2CRAM 37d5
	mkdir -p ${out_dir}
	ref=$4
	case $ref in
		37)
		# REF=/lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta
		REF=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
		;;
		37d5)
		# REF=/lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa
		REF=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
		;;
		38)
		# REF=/lustre/scratch114/resources/ref/Homo_sapiens/GRCh38_15/Homo_sapiens.GRCh38_15.fa #this should give an error
		# REF=/lustre/scratch113/projects/carl_seq/04072015/ref_check/Homo_sapiens.GRCh38_15.fa #this should give an error
		REF=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa #this should give an error
		;;
	esac
	filename=`basename ${file}`
	out_name="${filename%.*}"
	samtools view -h -T ${REF} ${file} -C -o ${out_dir}/${out_name}.cram
	
	#compare stats from bam and cram to check they're the same
	
	;;
	INDEX )
	# generate/overwrite bam index file
	# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/carl_seq/NEWBATCH/new_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_idx[1-${size}]" -o "LOGS/%J_bam_idx.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh /lustre/scratch113/projects/carl_seq/NEWBATCH/new_batch_files.list /lustre/scratch113/projects/carl_seq/NEWBATCH/ INDEX 
	# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/fvg_seq/NEWBATCH/new_batch_files.list|cut -f 1 -d " "`;bsub -J "bam_idx[1-${size}]" -o "LOGS/%J_bam_idx.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh /lustre/scratch113/projects/fvg_seq/NEWBATCH/new_batch_files.list /lustre/scratch113/projects/fvg_seq/NEWBATCH/ INDEX 
	mkdir -p ${out_dir}
	samtools index -b ${file}
	;;
	CHECKREF)
	#call GATK with a depth check command to check if there is an error with the reference selected
	module add hgi/gatk-protected/latest
	mkdir -p ${out_dir}
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
	mkdir -p ${out_dir}
	# use one of the sample bams present here, for VBI cohort: ~/esgi_vbseq/release-vb_sanger/20130807/sample_improved_bams_hgi_2/
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
	mkdir -p $out_dir/BAMCHECK_STATS/PLOTS/QC_GRIND/HTML

	echo "Processing $out_name"
	#create bam stats if they're not present already
	# /software/vertres/codebase/scripts/bamcheck -c 1,50,1 -d $1 > $out_dir/BAMCHECK_STATS/$out_name.bamchek.stats
	if [[ ! -s $out_dir/BAMCHECK_STATS/$out_name.bamchek.stats ]]
	then
		/software/vertres/codebase/scripts/bamcheck -d $1 > $out_dir/BAMCHECK_STATS/$out_name.bamchek.stats
	fi

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
	/software/hgi/pkglocal/samtools-0.1.19/bin/plot-bamcheck -p $out_dir/BAMCHECK_STATS/PLOTS/QC_GRIND/${out_name%%.*}/${out_name%%.*} $out_dir/BAMCHECK_STATS/${out_name}.bamchek.stats

	# generate a html report with images
	/nfs/users/nfs_m/mc14/Work/bash_scripts/bam_check_html_report.sh ${out_name%%.*} ../${out_name%%.*}/${out_name%%.*}-acgt-cycles.png ../${out_name%%.*}/${out_name%%.*}-coverage.png ../${out_name%%.*}/${out_name%%.*}-gc-content.png ../${out_name%%.*}/${out_name%%.*}-gc-depth.png ../${out_name%%.*}/${out_name%%.*}-indel-cycles.png ../${out_name%%.*}/${out_name%%.*}-indel-dist.png ../${out_name%%.*}/${out_name%%.*}-insert-size.png ../${out_name%%.*}/${out_name%%.*}-quals.png ../${out_name%%.*}/${out_name%%.*}-quals2.png ../${out_name%%.*}/${out_name%%.*}-quals3.png ../${out_name%%.*}/${out_name%%.*}-quals-hm.png > $out_dir/BAMCHECK_STATS/PLOTS/QC_GRIND/HTML/${out_name%%.*}_summary.html

  	;;
  	SEXCHECK)
	#Check sex in ba files
	mkdir -p ${out_dir}
	y_reads=`samtools view ${file} Y| cut -f 1 | sort -u | wc -l`
	echo "${file} ${y_reads}"
	;;

	FOFN )
	# mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/carl_seq/NEWBATCH/new_batch_files.list|cut -f 1 -d " "`;bsub -J "fofn[1-${size}]" -o "LOGS/%J_fofn.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/bam_check_pipeline.sh /lustre/scratch113/projects/carl_seq/NEWBATCH/new_batch_files.list /lustre/scratch113/projects/carl_seq/04072015/fofn FOFN carl_seq
	# generate a File of File names for bam improvement with the following header/info:
	# path lane center_name platform reads bases paired library sample study

	# file=603.rmdup.bam
	mkdir -p ${out_dir}
	REF=/lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa
	# /software/jre1.7.0_25/bin/java -Xmx2000m -Xms2000m -server -XX:+UseSerialGC -jar /software/hgi/pkglocal/picard-tools-1.127/picard.jar CollectAlignmentSummaryMetrics INPUT=${file} OUTPUT=${file}.stats STOP_AFTER=1000000 METRIC_ACCUMULATION_LEVEL=ALL_READS METRIC_ACCUMULATION_LEVEL=SAMPLE METRIC_ACCUMULATION_LEVEL=LIBRARY METRIC_ACCUMULATION_LEVEL=READ_GROUP 
	# /software/jre1.7.0_25/bin/java -Xmx2000m -Xms2000m -server -XX:+UseSerialGC -jar /software/hgi/pkglocal/picard-tools-1.127/picard.jar CollectAlignmentSummaryMetrics INPUT=${file} OUTPUT=${file}.align.stats STOP_AFTER=1000000 METRIC_ACCUMULATION_LEVEL=READ_GROUP
	# /software/jre1.7.0_25/bin/java -Xmx2000m -Xms2000m -server -XX:+UseSerialGC -jar /software/hgi/pkglocal/picard-tools-1.127/picard.jar CollectWgsMetrics INPUT=${file} OUTPUT=${file}.wgs.stats REFERENCE_SEQUENCE=${REF} STOP_AFTER=1000000
	
	# generate reads stats
	# samtools flagstat ${file} > ${file}.flagstat.stats
	path=${file}
	filename=`basename ${file}`
	lane=${filename%%.*}
	center_name=BGI
	platform=`samtools view -H ${file}|grep ^@RG|tr "\t" "\n"|grep PL|sed 's/PL://'|head -1`
	# reads=`samtools view -c ${file}`
	reads=`samtools idxstats ${file} | awk '{mapped +=$3;unmapped += $4}END{print mapped+unmapped}'`
	seq_length=`samtools view ${file}|head -1| awk '{print length($10)}'`
	bases=$[seq_length * reads]
	paired=YES
	library=`samtools view -H ${file}|grep ^@RG|tr "\t" "\n"|grep LB|sed 's/LB://'|head -1`
	sample=`samtools view -H ${file}|grep ^@RG|tr "\t" "\n"|grep SM|sed 's/SM://'|head -1`
	study=$4
	
	(echo "path lane center_name platform reads bases paired library sample study";echo "${file} ${lane} ${center_name} ${platform} ${reads} ${bases} ${paired} ${library} ${sample} ${study}" ) > ${out_dir}/${filename}.fofn

	;;

	IMPROVE)
	#piece of pipeline to perform bam improvement in a sanger-hgi style
	# we need a 4 step procedure: 
	# 0) Mark duplicates
	# 	- Use picard to mark, not remove, duplicates reads
	# 1) Realignement around INDELs
	# 	- Use GATK to perform local realignment of reads to correct misalignments due to the presence of indels. It's a two step procedure.
	# 		- First you define intervals for the Local Indel Realigner to target for realignment with the RealignerTargetCreator.
	# 		- Second you apply IndelRealigner walker to your data using the intervals previously generated

	# 2) Base Quality Recalibration
	# 	- Use GATK to perform a two-pass procedure to recalculate the base quality score
	# 		- Generates recalibration table with BaseRecalibrator based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).
	# 		- apply recalibration table
	
	# 3) MD tag recalculation
	# 	- Use samtools calmd to recarculate BAM's MD tag

	;;
esac
