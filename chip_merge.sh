#!/usr/local/bin/bash

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
awk 'FNR==NR { a[$1,$4,$5,$6]=$0; next } ($1,$4,$5,$6) in a { print $0,a[$1,$4,$5,$6] }' /netapp/dati/WGS_REF_PANEL/genotypes/VBI/shapeit/all_chr_VBI.bim /netapp/dati/WGS_REF_PANEL/genotypes/VBI/exome/M00705_plate1-19_1719_genotypes_zCall_converted.flipped.bim| tr " " "\t"| awk '{print $8,$2}' > exome_id_to_update.txt

for i in {1..22}
do
	mkdir -p /netapp/dati/WGS_REF_PANEL/genotypes/VBI/exome/${i}
	plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/VBI/exome/M00705_plate1-19_1719_genotypes_zCall_converted.flipped --update-name /netapp/dati/WGS_REF_PANEL/genotypes/VBI/exome/exome_id_to_update.txt 1 2 --chr ${i} --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/VBI/exome/${i}/chr${i}
done


#first create recoded files
for i in {1..22}
do
	mkdir -p /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/
	plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/chr${i} --recode --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}
	cut -f 2- -d " " /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}.ped | awk '{print $1,$0}' | sort -g -k1,1 > /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}_sorted.ped
	mv /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}.map /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}_sorted.map
	plink --file /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}_sorted --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/shapeit/${i}/chr${i}_sorted
done

for i in {1..22}
do
	mkdir -p /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/
	plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i} --recode --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}
	cut -f 2- -d " " /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}.ped | awk '{print $1,$0}' | sort -g -k1,1 > /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_sorted.ped
	mv /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}.map /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_sorted.map
	plink --file /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_sorted --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/exome/${i}/chr${i}_sorted
done
# then duplicate FID and IID
# and last, sort based on FID/IID

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

#clean genotypes removing samples with only exome data (low call rate)
pop=CARL
for i in {1..22}
do
    mkdir -p /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}
    plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/${i}/chr${i} --mind 0.1 --make-bed --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}
done

# calculate freqs for merged dataset
for pop in FVG CARL
do
for i in {1..22}
do
    plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i} --freq --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_freq
    awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_freq.frq > /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${i}_freq.tab
done
done
