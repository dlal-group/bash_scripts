#!/usr/local/bin/bash
#script to extract additional information for gemma results for UK10K replication
#Args:
#$1=gemma_pheno_list.txt path 
#$2=results path
#$3=imputed path files
#$4=sample file path in impute format (/nfs/users/nfs_m/mc14/Work/SANGER/FVG/PHENO/FVG.custom_sample)
#$5=out_path

if [ $# -lt 4 ]
then
    echo "missing arguments!!"
    echo "USAGE: gemma_additional_stats_pipeline.sh <gemma pheno list file path> <gemma results path> <imputed files path > <sample file>"
    echo "sample bsub command (ie. Lipids in FVG cohort):"
    echo "bsub -J \"additional_stats\" -o \"%J_additional_stats_LIPIDS.o\" -M4000 -R\"select[mem>4000] rusage[mem=4000]\" -q normal  -- gemma_additional_stats_pipeline.sh /nfs/users/nfs_m/mc14/Work/SANGER/FVG/PHENO/LIPIDS/new/gemma_pheno_list.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/LIPIDS/NEW_ORDER/output /lustre/scratch113/projects/uk10k/users/jh21/imputed/fvg/uk10k1kg.shapeit /nfs/users/nfs_m/mc14/Work/SANGER/FVG/PHENO/FVG.custom_sample" 
    exit 1
fi

#1) Extract for each trait the list of samples analized

out_dir=`dirname $1`;
res_path=$2;

while read line
do
    echo ${line}
    col=`echo ${line}| cut -f 1 -d " "`
    trait=`echo ${line}| cut -f 2 -d " "`
    
    echo "Col n: ${col}"
    echo "Trait: ${trait}"

    if [ ! -f ${out_dir}/${trait}.gemma_excludeindlist.txt ]; then
        awk -v col=${col} '$col=="NA"' ${out_dir}/gemma_pheno.txt | cut -f 1 -d " " >  ${out_dir}/${trait}.gemma_excludeindlist.txt
        nind=`wc -l ${out_dir}/${trait}.gemma_excludeindlist.txt | cut -f 1 -d " "`
        echo "N analized ind: ${nind}"
    fi
    
    for chr in {1..22} X
    # for chr in X
    do
    	echo "Chr ${chr}"
	    if [ ! -s ${res_path}/${trait}.${chr}.gemma_snplist.txt ]; then
	            #Extract from the result file the list of snp you want to include in the analysis
	            #check what type of result files we are looking at
	        if [ $chr == "X" ]
	        then
	        	chr_code=`cut -f 1 ${res_path}/${trait}.chrX.tab.assoc.txt | tail -n+2 | sort | uniq`

				if [[ "$chr_code" != "X" ]]
				then				
					sed 's/^'$chr_code'	/23	/g' ${res_path}/${trait}.chrX.tab.assoc.txt | gzip -c > ${res_path}/${trait}.chrX.tab.assoc.txt.gz
				fi
	            zgrep "^${chr}	" ${res_path}/${trait}.chrX.tab.assoc.txt.gz|cut -f 3 > ${res_path}/${trait}.${chr}.gemma_snplist.txt
	        else
	            grep "^${chr}	" ${res_path}/${trait}.chr${chr}.tab.assoc.txt |cut -f 3 > ${res_path}/${trait}.${chr}.gemma_snplist.txt
	        fi
	    else
	    	echo "File ${res_path}/${trait}.${chr}.gemma_snplist.txt already existing!!"
	    fi
       gemma_additional_stats_calc_snptest_f3_fvg.sh $3 $4 ${res_path}/${trait}.${chr}.gemma_snplist.txt ${out_dir}/${trait}.gemma_excludeindlist.txt ${res_path}/added_stats ${chr}
    done

done < <(tail -n+2 $1)
