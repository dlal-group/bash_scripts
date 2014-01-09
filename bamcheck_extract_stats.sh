#!/usr/local/bin/bash

#extract info per sample from bamcheck results
#Args
#$1 = stats files path

# files=(`ls $1/*.bamchek.stats`)
# index=$[LSB_JOBINDEX]

# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
out_name=`basename $1`

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
plot-bamcheck -r /lustre/scratch111/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa.gc_stats -p $2/PLOTS/QC_GRIND/${out_name%%.*}/ $1
