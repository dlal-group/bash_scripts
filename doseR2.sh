#fast script to calculate R squared correlation between genotypes and imputed genotypes
#args
chr=$1
wgs_gen_file_path=$2
imp_gen_file_path=$3
output_path=$4

#create outdir
mkdir -p ${output_path}

#count how many samples (all files must have the same sample number)
n_sam1=`zcat ${wgs_gen_file_path}| awk '{print (NF-5)/3}'| sort| uniq`
n_sam2=`zcat ${imp_gen_file_path}| awk '{print (NF-5)/3}' |sort|uniq`

echo $n_sam1
echo $n_sam2

if [ $n_sam1 -eq $n_sam2 ]
then
	n_sam=$n_sam1
else
	echo "ERROR!! the files must have the same sample size!!"
	exit 1
fi

#convert gen dosages to mach format dosages for the imput files
# zcat ${wgs_gen_file_path} | awk -v chr=${chr} -v n_sam=${n_sam} '{printf "chr"chr":"$3; for(i=1;i<=n_sam;i++) printf " " $(i*3+3)*2+$(i*3+4); printf "\n" }' > ${output_path}/chr${chr}.wgs.dose
zcat ${wgs_gen_file_path} | awk -v chr=${chr} -v n_sam=${n_sam} '{printf "chr"chr":"$3" "$4" "$5; for(i=1;i<=n_sam;i++) printf " " $(i*3+3)*2+$(i*3+4); printf "\n" }' > ${output_path}/chr${chr}.wgs.dose

# zcat ${imp_gen_file_path} | awk -v chr=${chr} -v n_sam=${n_sam} '{printf "chr"chr":"$3; for(i=1;i<=n_sam;i++) printf " " $(i*3+3)*2+$(i*3+4); printf "\n" }' > ${output_path}/chr${chr}.imp.dose
zcat ${imp_gen_file_path} | awk -v chr=${chr} -v n_sam=${n_sam} '{printf "chr"chr":"$3" "$4" "$5; for(i=1;i<=n_sam;i++) printf " " $(i*3+3)*2+$(i*3+4); printf "\n" }' > ${output_path}/chr${chr}.imp.dose
 
#convert tped plink format file 
#cat ${ped_file} | sed 's/\t1 1/\t0/g' | sed 's/\t1 2/\t1/g' | sed 's/\t2 2/\t2/g' > ${output_path}/chr${chr}.wgs.dose

#find matches between files
#awk '{ if (FILENAME ~/wgs/) {$1=$3=$4=""; pos[$2]=$0} else { if (pos[$1] !="") print pos[$1],$0 } }' ${output_path}/chr${chr}.wgs.dose ${output_path}/chr${chr}.imp.dose > ${output_path}/chr${chr}.matched
awk '{ if (FILENAME ~/wgs/) {pos[$1]=$0} else { if (pos[$1] !="") print pos[$1],$0 }}' ${output_path}/chr${chr}.wgs.dose ${output_path}/chr${chr}.imp.dose > ${output_path}/chr${chr}.matched

#corrected the eventually flipped alleles and skip totally mismatching alleles
# awk -v n_sam=${n_sam} '{
# if($2 == $(n_sam+5) && $3 == $(n_sam+6))
#   {printf $1;for(i=4;i<=(n_sam+3);i++) printf " " $i;printf " " $(n_sam+4);for(j=(n_sam+7);j<=(2*(n_sam+3));j++)printf " "$j;printf "\n";}
# else if ($3 == $(n_sam+5) && $2 == $(n_sam+6))
#   {printf $1;for(i=4;i<=(n_sam+3);i++) printf " " 2-$i;printf " " $(n_sam+4);for(j=(n_sam+7);j<=(2*(n_sam+3));j++)printf " "$j;printf "\n";}
# }' ${output_path}/chr${chr}.matched > ${output_path}/chr${chr}.matched_clean

awk -v n_sam=110 '{
if($2 == $(n_sam+5) && $3 == $(n_sam+6))
  {printf $1;for(i=4;i<=(n_sam+3);i++) printf " " $i;printf " " $(n_sam+4);for(j=(n_sam+7);j<=(2*(n_sam+3));j++)printf " "$j;printf "\n";}
}' chr7.matched > chr7.test

#use matching sites to calculate the R2 correlation
cat ${output_path}/chr${chr}.matched_clean | awk -v n_sam=${n_sam} '{ sx=sy=sxy=sx2=sy2=0; for (i=2;i<=(n_sam+1);i++){ sx+=$i; sy+=$(i+n_sam+1); sxy+=$i*$(i+n_sam+1); sx2+=$i*$i; sy2+=$(i+n_sam+1)*$(i+n_sam+1) } if ( (n_sam*sx2-sx*sx)*(n_sam*sy2-sy*sy) !=0 ) print $1, (n_sam*sxy-sx*sy)*(n_sam*sxy-sx*sy)/((n_sam*sx2-sx*sx)*(n_sam*sy2-sy*sy));else print $1,"NA" }' > ${output_path}/chr${chr}.doseR2
# cat ${output_path}/chr${chr}.matched_clean | awk -v n_sam=${n_sam} '{ sx=sy=sxy=sx2=sy2=0; for (i=2;i<=(n_sam+1);i++){ sx+=$i; sy+=$i; sxy+=$i*$i; sx2+=$i*$i; sy2+=$i*$i } if ( (n_sam*sx2-sx*sx)*(n_sam*sy2-sy*sy) !=0 ) print $1, (n_sam*sxy-sx*sy)*(n_sam*sxy-sx*sy)/((n_sam*sx2-sx*sx)*(n_sam*sy2-sy*sy));else print $1,"NA" }' > ${output_path}/chr${chr}.doseR2

