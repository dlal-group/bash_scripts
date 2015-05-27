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
