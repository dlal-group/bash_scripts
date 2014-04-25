#!/usr/local/bin/bash
##############ATTENTION!!!!!################################
#This script works with COMPRESSED FILES in PLINK FORMAT!!##
############################################################
#We can use this same launcher to check data using BASH or R
#We will pass a flag for the R or BASH way
#Args:
#$1 = file to process
#$2 = gz flag
#$3 = trait
if [ $# -lt 1 ]
then
	echo -e "##############ATTENTION!!!!!################################\n#This script works with COMPRESSED FILES in PLINK FORMAT!!##\n############################################################"
	echo "MISSING ARGUMENTS!!!"
	echo -e "USAGE:"
	echo -e "BASH mode:\ncheck_pvals.sh <file_path>"
	echo -e "R mode:\ncheck_pvals.sh <file_path> <gzip flag> <trait name>"
	echo -e "<file_path>: file with results to check\n<gzip flag>: [TRUE/FALSE] used to tell R if the file is in gz format\n<trait name> : name of the trait"
	exit 1
fi


if [ $# -lt 2 ]
then
	echo "Shell Mode"
	echo -e "Processing file $1"
	#check pval calculating lambda and qqplot
	# (echo "chr	rs	pos	chi2	pval";cut -f -3,5,8,9 $1 | tail -n+2 | awk '{OFS=" "}{print $0,(($4/$5)*($4/$5))}' ) | tr " " "\t" | gzip -c > $1_chi2.gz
	# (echo "chr	rs	pos	chi2	pval";zcat $1| sed 's/ \+/ /g' | sed 's/^ //g' | tr " " "\t" |cut -f -3,8,9 | tail -n+2 )| tr " " "\t" | gzip -c > $1_chi2.gz
	(echo "chr	rs	pos	chi2	pval";zcat $1 |cut -f -3,8,9 | tail -n+2 )| tr " " "\t" | awk '{if($5 != "NA") print $0}'| gzip -c > $1_chi2.gz

	#we can make qqplots with gnuplots!
	#create the file to plot a qqplot and the manhattan plot
	#count rows of the file
	t=`zcat $1_chi2.gz | tail -n+2|wc -l`
	(echo "chr rs pos chi2 pval -log_P exp_p -log_exp_p"|tr " " "\t";(zcat $1_chi2.gz | tail -n+2 | awk '{print $0,-(log($5)/log(10.0))}'|tr " " "\t"|sort -g -k5,5| awk -v l=$t '{m=l;a=0.5;print $0,((NR)-a)/(m+(1-a)-a),-(log(((NR)-a)/(m+(1-a)-a))/log(10.0))}'|tr " " "\t"))> $1_chi2.plots
	# lets try to calculate the same formula with awk : lambda=median(chi2,na.rm=T)/qchisq(0.5,1)
	# qchisq(0.5,1) = 0.4549364 
	#This calculate the MEDIAN value in a NUMERICALLY SORTED column without NA:
	# awk '{count[NR] = $2;}END{if (NR % 2) {print count[(NR + 1) / 2];} else {print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0;}}' 

	md=`zcat $1_chi2.gz | tail -n+2 | cut -f 4 | fgrep -v "NA" | sort -g | awk '{count[NR] = $1}END{if (NR % 2) {print count[(NR + 1) / 2];} else {print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0;}}'`
	median=`printf "%f\n" $md`
	l=$( bc -l <<< "${median}/0.4549364")
	lambda=$(echo "$l")
	check=$( bc -l <<< $l'>'1.01)
	echo "#######################################"
	echo "The inflation factor fo this trait is:"
	echo "$l"

	if [ $check -eq 1 ]
	then
		echo "You need to correct your p-values!!"
		#calculate the corrected chi2
		(echo "rs chi2_gc";zcat $1_chi2.gz | tail -n+2 | awk -v lam=$lambda '{print $2, $4/lam}')| tr " " "\t" | gzip -c > $1_chi2_gc.gz
		echo "Corrected chi^2 in $1\_chi2_gc.gz"
	else
		echo "You don't need to correct your p-values!!"
	fi
	echo "#######################################"
else
	echo "R mode!"
	if [ $# -lt 3 ]
	then
		echo "MISSING ARGUMENTS for R mode!!!"
		echo -e "USAGE:"
		echo -e "R mode:\ncheck_pvals.sh <file_path> <gzip flag> <trait name>"
		echo -e "<file_path>: file with results to check\n<gzip flag>: [TRUE/FALSE] used to tell R if the file is in gz format\n<trait name> : name of the trait"
		exit 1
	fi
	echo -e "Processing file $1\non trait $3"

	#check pval calculating lambda and qqplot
	(echo "chr	rs	pos	beta	se	p_wald	chi2";cut -f -3,5,6,9 $1 | tail -n+2 | awk '{OFS=" "}{print $0,(($4/$5)*($4/$5))}' ) | tr " " "\t"| gzip -c > $1_chi2.gz
	#Call the R script
	R CMD BATCH '--args '$1_chi2.gz' 'TRUE' '$3'' ~/Work/r_scripts/result_check.r $3.Rout
fi

