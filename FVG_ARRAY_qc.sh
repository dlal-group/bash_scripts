#####snippets to work on genotype qc
#Save each output step in a corrispondent folder to avoid caos
# 1) Remove bad samples using PLINK:
# - convert report to ped, using convertReportToTPED.py

zcall_path=$1
raw_data_path=$2
gs_report=$3
prefix=`date +"%d%m%Y%H%M%S"`
pop=$4


zcall_path="/netapp/nfs/softwares/zCall/Version3_GenomeStudio/GenomeStudio"
raw_data_path="/home/shared/14022017_FVG_MEGA_QC/fvg_unico_mega_V1_2016-05-05"
gs_report=${raw_data_path}/GIANT_V1_clean.txt
gs_report_basename=`basename ${gs_report}`
prefix=`date +"%d%m%Y%H%M%S"`
pop="FVG"

mkdir -p ${raw_data_path}/01_het_check
mkdir -p ${raw_data_path}/02_callrate_check
mkdir -p ${raw_data_path}/03_cluster_check


${zcall_path}/convertReportToTPED.py -R ${gs_report} -O ${raw_data_path}/01_het_check/${prefix}_${pop}

# check het and missingness with plink

plink --tfile ${raw_data_path}/01_het_check/${prefix}_${pop} --het --out ${raw_data_path}/01_het_check/${prefix}_${pop}_het_check

# add het column
tail -n+2 ${raw_data_path}/01_het_check/${prefix}_${pop}_het_check.het| awk '{OFS="\t"}{print $0,($5-$3)/$5}'|sort -g -r -k7,7 > ${raw_data_path}/01_het_check/${prefix}_${pop}_het_check.het.rate

#call rate check
# extract samples with less tha 0.95 call rate
plink --tfile ${raw_data_path}/01_het_check/${prefix}_${pop} --mind 0.05 --make-bed --out ${raw_data_path}/02_callrate_check/${prefix}_${pop}_mind005_check

# check if those samples are the same with high het-rate
fgrep -w -f <(cut -f 1 ${raw_data_path}/02_callrate_check/${prefix}_${pop}_mind005_check.irem) ${raw_data_path}/01_het_check/${prefix}_${pop}_het_check.het.rate > ${raw_data_path}/02_callrate_check/HET_callrate_bad_samples.list

# remove bad samples from report
/netapp/nfs/softwares/zCall/additionalScripts/dropSamplesFromReport_FasterVersion.py ${gs_report} ${raw_data_path}/02_callrate_check/${prefix}_${pop}_mind005_check.irem > ${raw_data_path}/02_callrate_check/${gs_report_basename}_callrate95.txt

#Use the script findMeanSD.py to calculate μ and σ of both homozygote clusters for common sites (MAF > 5%)
#we need to use the cleaned report
${zcall_path}/findMeanSD.py -R ${gs_report} > ${raw_data_path}/03_cluster_check/${gs_report_basename}_mean_sd.txt


Take the output of Step 1 and run findBetas.r to derive the linear regression model
	We suggest using the “1” flag for weighted linear regression

Use the output from findBetas.r and the clean GSR to determine which z-score threshold to use by running calibrateZ.py for different z-score threshold inputs (i.e. 3-15)
There was a small bug in the calibrateZ.py script making it so the script would not run.

Determine which z-score threshold works best for your dataset based on the output statistics

Use the findThresholds.py script to derive the thresholds using the linear regression model, the optimal z-score threshold, and the clean GSR


Use the zCall.py script with the clean (or original) GSR and the thresholds output from findThresholds.py. Output is a TPED and TFAM file with only No Calls recalled by zCall.
