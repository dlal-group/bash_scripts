#####snippets to work on genotype qc

# 1) Remove bad samples using PLINK:
# - convert report to ped, using convertReportToTPED.py

zcall_path="/netapp/nfs/softwares/zCall/Version3_GenomeStudio/GenomeStudio"

${zcall_path}/convertReportToTPED.py -R GIANT_V1.txt -O /home/shared/14022017_FVG_MEGA_QC/fvg_unico_mega_V1_2016-05-05/14022017_FVG_MEGA_V1

# check het and missingness with plink

plink --tfile /home/shared/14022017_FVG_MEGA_QC/fvg_unico_mega_V1_2016-05-05/14022017_FVG_MEGA_V1 --het --out /home/shared/14022017_FVG_MEGA_QC/fvg_unico_mega_V1_2016-05-05/14022017_FVG_MEGA_V1_het_check

# add het column
tail -n+2 14022017_FVG_MEGA_V1_het_check.het| awk '{OFS="\t"}{print $0,($5-$3)/$5}' > 14022017_FVG_MEGA_V1_het_check.het.rate

# extract samples with less tha 0.95 call rate
plink --tfile /home/shared/14022017_FVG_MEGA_QC/fvg_unico_mega_V1_2016-05-05/14022017_FVG_MEGA_V1 --mind 0.05 --make-bed --out /home/shared/14022017_FVG_MEGA_QC/fvg_unico_mega_V1_2016-05-05/14022017_FVG_MEGA_V1_mind001_check

# check if those samples are the same with high het-rate
fgrep -w -f <(cut -f 1 14022017_FVG_MEGA_V1_mind001_check.irem) 14022017_FVG_MEGA_V1_het_check.het.rate

# remove bad samples from report
/netapp/nfs/softwares/zCall/additionalScripts/dropSamplesFromReport_FasterVersion.py -R GIANT_V1.txt

Use the script findMeanSD.py to calculate μ and σ of both homozygote clusters for common sites (MAF > 5%)

Take the output of Step 1 and run findBetas.r to derive the linear regression model
	We suggest using the “1” flag for weighted linear regression

Use the output from findBetas.r and the clean GSR to determine which z-score threshold to use by running calibrateZ.py for different z-score threshold inputs (i.e. 3-15)
There was a small bug in the calibrateZ.py script making it so the script would not run.

Determine which z-score threshold works best for your dataset based on the output statistics

Use the findThresholds.py script to derive the thresholds using the linear regression model, the optimal z-score threshold, and the clean GSR


Use the zCall.py script with the clean (or original) GSR and the thresholds output from findThresholds.py. Output is a TPED and TFAM file with only No Calls recalled by zCall.
