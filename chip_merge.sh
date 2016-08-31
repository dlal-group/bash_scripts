#!/usr/bin/env bash

#script to merge exome chip and array chip together

pop=VBI
#prepare FVG and VBI genotype files
for i in {1..22}
do
mkdir -p ${i}
gunzip chr${i}.hap.gz 
# mv chr${i}.hap 
cut -f 2- -d " " chr${i}.hap | awk -v chr=${i} '{print chr,$0}' > chr${i}.haps
gzip -c chr${i}.haps > chr${i}.hap.gz
/netapp/nfs/Max/software/shapeit/bin/shapeit -convert --input-haps chr${i} --output-vcf chr${i}.vcf
plink --vcf chr${i}.vcf --chr ${i} --make-bed --out ${i}/chr${i}
rm chr${i}.vcf
done

#map exome chip sites to array variants
# awk 'FNR==NR { a[$1,$4,$5,$6]=$0; next } ($1,$4,$5,$6) in a { print $0,a[$1,$4,$5,$6] }' /netapp/dati/WGS_REF_PANEL/genotypes/VBI/shapeit/all_chr_VBI.bim /netapp/dati/WGS_REF_PANEL/genotypes/VBI/exome/M00705_plate1-19_1719_genotypes_zCall_converted.flipped.bim| tr " " "\t"| awk '{print $8,$2}' > exome_id_to_update.txt

for i in {1..22}
do
	mkdir -p /netapp/dati/WGS_REF_PANEL/genotypes/VBI/exome/${i}
	plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/VBI/exome/VBI_exome_name_updated_filtered_cleaned --chr ${i} --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/VBI/exome/${i}/chr${i}
done 


#first create recoded files: ARRAY
pop=VBI
pop=CARL
pop=FVG
for i in {1..22}
do
	mkdir -p /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/
	plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i} --recode --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}
	cut -f 2- -d " " /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}.ped | awk '{print $1,$0}' | sort -g -k1,1 > /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}_sorted.ped
	mv /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}.map /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}_sorted.map
	plink --file /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}_sorted --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}_sorted
done

#EXOME CHIP
pop=VBI
pop=CARL
pop=FVG
for i in {1..22}
do
	mkdir -p /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/
	plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i} --recode --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}
	cut -f 2- -d " " /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}.ped | awk '{print $1,$0}' | sort -g -k1,1 > /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_sorted.ped
	sed 's/exm-rs/rs/g' /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}.map > /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_sorted.map
	plink --file /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_sorted --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_sorted
done
# then duplicate FID and IID
# and last, sort based on FID/IID

pop=VBI
pop=CARL
pop=FVG

#merge data
# for pop in CARL FVG VBI
for pop in VBI
do
for i in {1..22}
do
    mkdir -p /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}
    if [ -s /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}/chr${i}-merge.missnp ]
    then
        plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_sorted --flip /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}/chr${i}-merge.missnp --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_flipped
        plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}_sorted --bmerge /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_flipped --merge-mode 2 --merge-equal-pos --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}/chr${i}
    else
        plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}_sorted --bmerge /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_sorted --merge-mode 2 --merge-equal-pos --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}/chr${i}
    fi
done
done

#take care of snp in chr14 for VBI cohort
for pop in VBI
do
for i in 14
do
    mkdir -p /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}
    if [ -s /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}/chr${i}-merge.missnp ]
    then
        plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_sorted --flip /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}/chr${i}-merge.missnp --exclude-snps rs7824 --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_flipped
        plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}_sorted --bmerge /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_flipped --merge-mode 2 --merge-equal-pos --exclude-snps rs7824 --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}/chr${i}
    else
        plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}_sorted --bmerge /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_sorted --merge-mode 2 --merge-equal-pos --exclude-snps rs7824 --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}/chr${i}
    fi
done
done

#clean genotypes removing problematic sites
for pop in CARL FVG VBI
for pop in VBI
do
fgrep "Multiple positions" /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/*/*.log|awk '{print $7}'| sed "s/'//g"|sed 's/\.//g' > /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/snps_to_remove.list
fgrep "same position" /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/*/*.log|awk '{print $3,$5}'|tr " " "\n"| sed "s/'//g" >> /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/snps_to_remove.list

for i in {1..22}
do
plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}/chr${i} --exclude /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/snps_to_remove.list --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}/chr${i}_cleaned
done
done


#clean genotypes removing samples with only exome data (low call rate)
for pop in CARL FVG VBI
do
for i in {1..22}
do
    mkdir -p /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}
    plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}/chr${i}_cleaned --mind 0.1 --hwe 0.00001 --maf 0.0009 --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}
done
done

#check merged genotypes
for pop in CARL FVG VBI
do
for i in {1..22}
do
head -1 /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}.bim
done
done


# calculate freqs for merged dataset
for pop in FVG CARL VBI
do
for i in {1..22}
do
    plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i} --nonfounders --freq --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_freq
    awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_freq.frq > /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_freq.tab
done
done

# print all freq together
for pop in FVG CARL VBI
do
(echo -e "CHR\tSNP\tA1\tA2\tMAF\tNCHROBS";
for i  in {1..22}
do
fgrep -v NCHROBS /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_freq.tab
done) > /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${pop}_freq.tab
done

#run phasing for all merged chromosomes
#submit jobs using pbs
for pop in FVG CARL VBI
do
for i in {1..22}
do
echo "/netapp/nfs/Max/software/shapeit/bin/shapeit -B /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i} -M /netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr${i}_combined_b37.txt -O /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}.haps.gz /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}.sample -T 8"| qsub -o /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_shapeit.log -e /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_shapeit.e -V -N ${pop}_chr${i}_shapeit -pe openmpi-fillup 1-8
done
done

####################################################################################
#11/08/2016
#Prepare input data for Impute test on MATULLO's samples
#clean genotypes removing samples with only exome data (low call rate)
for pop in MATULLO
do
for i in 2
do
    mkdir -p /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}
    mkdir -p /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}
    plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/${i}/chr${i} --mind 0.01 --hwe 0.00001 --geno 0.01 --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}/chr${i}
    #clean all duplicate position
    plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}/chr${i} --exclude /netapp/dati/WGS_REF_PANEL/genotypes/MATULLO/${i}/chr${i}_dupe_rs.list --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}
done
done

# calculate freqs for merged dataset
for pop in MATULLO
do
for i in 2
do
    plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i} --nonfounders --freq --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_freq
    awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_freq.frq > /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_freq.tab
done
done

# print all freq together
for pop in MATULLO
do
(echo -e "CHR\tSNP\tA1\tA2\tMAF\tNCHROBS";
for i  in 2
do
fgrep -v NCHROBS /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_freq.tab
done) > /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/${pop}_freq.tab
done

#run phasing for all merged chromosomes
for pop in MATULLO
do
for i in 2
do
echo "/home/cocca/softwares/bin/shapeit -B /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i} -M /netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr${i}_combined_b37.txt -O /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}.haps.gz /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}.sample -T 8"| qsub -o /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_shapeit.log -e /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_shapeit.e -V -N ${pop}_chr${i}_shapeit -l h_vmem=8G 
done
done


