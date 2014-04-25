#!/usr/local/bin/bash

#extract positive controls from results

#merge chunk files and harmonize snp rs id to be the same as bimbam format:
if [[ ${LSB_JOBINDEX} == "23" ]]
then
	c="X"
else
	c=${LSB_JOBINDEX}
fi

poscon=$1
res_path=$2
trait=$3
out=$4

#extract poscon for the trait category for each chr

# (head -1 ${poscon};awk -v chr=${c} '$4==chr' ${poscon}|sort -g -k5,5 | fgrep "${trait}") > ${out}/${trait}.${c}.poscon
if [[ ${trait} == "TC" || ${trait} == "LDL" || ${trait} == "HDL" || ${trait} == "TG" ]]
then
	echo "${trait} - category: LIPIDS"
	(head -1 ${poscon};awk -v chr=${c} '{if($4==chr && ($6~"TC" || $6~"LDL" || $6~"HDL" || $6~"TG")) print $0}' ${poscon}|sort -g -k5,5 ) > ${out}/${trait}.${c}.poscon
fi
if  [[ ${trait} == "GLU" || ${trait} == "GLUadjBMI" || ${trait} == "INS" || ${trait} == "INSadjBMI" || ${trait} == "HOMA_B" || ${trait} == "HOMA_BadjBMI" || ${trait} == "HOMA_IR" || ${trait} == "HOMA_IRadjBMI" ]]
then
	echo "${trait} - category: GLYCEMIC"
	(head -1 ${poscon};awk -v chr=${c} '{if($4==chr && ($6~"Glucose" || $6~"Insulin" || $6~"InsulinRes" || $6~"Insulin_BMIadj" || $6~"InsulinRes" || $6~"HOMA_b" || $6~"HOMA_b_BMIadj" || $6~"HOMA_ir" || $6~"HOMA_ir_BMIadj")) print $0}' ${poscon}|sort -g -k5,5 ) > ${out}/${trait}.${c}.poscon
fi

if [[ ${trait} == "HB" || ${trait} == "MCH" || ${trait} == "MCHC" || ${trait} == "MCV" || ${trait} == "PCV" || ${trait} == "PLT" || ${trait} == "RBC" || ${trait} == "WBC" ]]
then
	echo "${trait} - category: BLOOD"
	(head -1 ${poscon};awk -v chr=${c} '{if($4==chr && ($6~"HGB" || $6~"MCH" || $6~"MCHC" || $6~"MCV" || $6~"PCV" || $6~"PLT" || $6~"RBC" || $6~"WBC" || $6~"Platelet" || $6~"MPV" || $6~"Coagulation" || $6~"Hematocrit" || $6~"Hematological" || $6~"Hemo")) print $0}' ${poscon}|sort -g -k5,5 ) > ${out}/${trait}.${c}.poscon
fi
if [[ ${trait} == "BMI" || ${trait} == "HEIGHT" || ${trait} == "HIP" || ${trait} == "HIPadjBMI" || ${trait} == "TFM" || ${trait} == "TLM" || ${trait} == "WAIST" || ${trait} == "WAISTadjBMI" || ${trait} == "WEIGHT" || ${trait} == "WHR" || ${trait} == "WHRadjBMI" ]]
then
	echo "${trait} - category: ANTHROP"
	(head -1 ${poscon};awk -v chr=${c} '{if($4==chr && ($6~"BMI" || $6~"Height" || $6~"Obesity" || $6~"TFM" || $6~"TLM" || $6~"Waist" || $6~"Weight" || $6~"WHR" )) print $0}' ${poscon}|sort -g -k5,5 ) > ${out}/${trait}.${c}.poscon
fi
if [[ ${trait} == "DBP" || ${trait} == "HR" || ${trait} == "SBP" ]]
then
	echo "${trait} - category: CARDIO"
	(head -1 ${poscon};awk -v chr=${c} '{if($4==chr && ($6~"DBP" || $6~"HR" || $6~"SBP" || $6~"heartratepulse" || $6~"RR" || $6~"HT" )) print $0}' ${poscon}|sort -g -k5,5 ) > ${out}/${trait}.${c}.poscon
fi

if [[ ${trait} == "CRP" ]]
then
	echo "${trait} - category: INFL"
	(head -1 ${poscon};awk -v chr=${c} '{if($4==chr && ($6~"CRP" || $6~"Inflammatory")) print $0}' ${poscon}|sort -g -k5,5 ) > ${out}/${trait}.${c}.poscon
fi

if [[ ${trait} == "ALK" || ${trait} == "BIL" || ${trait} == "GGT" ]]
then
	echo "${trait} - category: LIVER"
	(head -1 ${poscon};awk -v chr=${c} '{if($4==chr && ($6~"Alkphosp" || $6~"Bilirubin" || $6~"GGT" || $6~"Albumin" )) print $0}' ${poscon}|sort -g -k5,5 ) > ${out}/${trait}.${c}.poscon
fi
if [[ ${trait} == "CREA" || ${trait} == "SOD" || ${trait} == "UREA" || ${trait} == "URICA" ]]
then
	echo "${trait} - category: RENAL"
	(head -1 ${poscon};awk -v chr=${c} '{if($4==chr && ($6~"Creatinine" || $6~"Urea" || $6~"UricAcid" || $6~"Renal" || $6~"Phosphate")) print $0}' ${poscon}|sort -g -k5,5 ) > ${out}/${trait}.${c}.poscon
fi


#extract for each chr results, the info for pos control
# (head -1 ${res_path}/${trait}.chr${c}.tab.assoc.txt;fgrep -w -f <( cut -f 5 -d " " ${out}/${trait}.${c}.poscon) ${res_path}/${trait}.chr${c}.tab.assoc.txt | sort -g -k3,3 ) | cut -f -6,9,11 > ${out}/${trait}.${c}.poscon.result
(head -1 ${res_path}/${trait}.chr${c}.tab.assoc.txt;fgrep -w -f <( cut -f 5 -d " " ${out}/${trait}.${c}.poscon) ${res_path}/${trait}.chr${c}.tab.assoc.txt | sort -g -k3,3 ) | cut -f -6,9,11 > ${out}/${trait}.${c}.poscon.result

#create a joint resume table with 
#CHR POS rsID P_repl beta se P_pubbl
(echo "CHR RS POS BETA SE P_wald P_score P_pval TRAIT";awk 'FNR==NR { a[$5]=$0; next } $3 in a { print $0,a[$3] }' ${out}/${trait}.${c}.poscon ${out}/${trait}.${c}.poscon.result |  tr "\t" " " | cut -f -3,5-8,11,14 -d " ") > ${out}/${trait}.${c}.poscon.join

