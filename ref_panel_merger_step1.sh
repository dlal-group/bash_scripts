#!/usr/local/bin/bash

#script to create ref panel from different populations suitables for merging with impute
#Args:
#$1= basepath first panel input files
#$2= basepath second panel input files
#$3= main output folder
#$4=POP 1 name
#$5=POP 2 name
#$6= indels files path


first_panel_path=$1
second_panel_path=$2
output_basepath=$3
pop1=$4
pop2=$5
indels_files_path=$6



function build_panel_input_template(){
cat << EOF
#!/usr/local/bin/bash
#build the panel input script
#we need as inputs:
#\$1=pop 1 (vbi) SNPs vcf filename =>$1
#\$2=pop 1 (vbi) INDELS vcf filename =>$2
#\$3=pop 2 vcf filename =>$3
#\$4=pop 1 output path =>$4
#\$5=pop 2 output path =>$5
#\$6=merged output =>$6
#\$7=chr
#\$8=pop 1
#\$9=pop 2

first_panel_file_path=$1
indels_file_path=$2
second_panel_file_path=$3
pop1_outpath=$4
pop2_outpath=$5
merged_outpath=$6
chr=$7
pop1=$8
pop2=$9

#set UK10K AC for singleton with REF allele
if [ \${pop2} == "UK10K" ]
then
	SAC=7561
fi

if [ \${pop2} == "1000GP" ]
then
	SAC=2183
fi

#merge back snps and indels
#we can keep all separated..'cause it's convenient for us
#vcf-concat /lustre/scratch110/sanger/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr7.re_ann.vcf.gz /lustre/scratch110/sanger/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/INDELS/CHR7/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.INDELS.7.vcf.gz | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.chr7.vcf.gz

#those first lines are only appliable to the VBI population because of the separation between indels and snps files
#convert indels vcf to tab format. Remove also those sites that are multiallelic in the vcf file
echo "Convert \${pop1} files in tab format..."
vcf-query \${indels_file_path} -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\n" |tr "\t" " " | awk '{if(\$3 == ".") print \$1,\$2,"chr"\$1":"\$2,\$4,\$5,\$6,\$7;else print \$0}' | awk '!(\$5~",")' > \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.tab

#convert SNPS files to tab format
vcf-query \${first_panel_file_path} -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\n" |tr "\t" " "| awk '{if(\$3 == ".") print \$1,\$2,"chr"\$1":"\$2,\$4,\$5,\$6,\$7;else print \$0}' | awk '!(\$5~",")' > \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.tab

echo "DONE!"

echo "Convert \${pop2} files in tab format..."
#same for 1000G panel
vcf-query \${second_panel_file_path} -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\n"| tr "\t" " " | awk '{if(\$3 == ".") print \$1,\$2,"chr"\$1":"\$2,\$4,\$5,\$6,\$7;else print \$0}' | awk '!(\$5~",")' > \${pop2_outpath}/\${pop2}.chr\${chr}.tab

echo "DONE!"
######(FOR SNPs ONLY)########
#once we have the tabbed files we check for allele mismatch in overlapping sites between VBI and other panels 
echo "Check for alleles mismatches in SNPS..."
#Same alleles
awk 'NR==FNR{a[\$2,\$4,\$5]=\$0;next;}(a[\$2,\$4,\$5] && length(\$4)==1 && length(\$5)==1)' \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.tab \${pop2_outpath}/\${pop2}.chr\${chr}.tab > \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching.keeplist 

#flipped alleles: we are almost sure that there are not flipped alleles!!...but we need to check!
awk 'NR==FNR{a[\$2,\$4,\$5]=\$0;next;}(a[\$2,\$5,\$4] && length(\$4)==1 && length(\$5)==1)' \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.tab \${pop2_outpath}/\${pop2}.chr\${chr}.tab > \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_flipped.keeplist
if [ -s \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_flipped.keeplist ]
then
  echo "##############################"
  echo "WARNING!!!FLIPPED GENOTYPES!!!"
  echo "WARNING!!!FLIPPED GENOTYPES!!!"
  echo "WARNING!!!FLIPPED GENOTYPES!!!"
  echo "##############################"
  exit 1
fi
#awk 'NR==FNR{a[\$2,\$4,\$5]=\$0;next;}(a[\$2,\$5,\$4] && length(\$4)==1 && length(\$5)==1)' \${first_panel_path}/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr7.re_ann.tab \${second_panel_path}/ALL.chr7.phase1_release_v3.20101123.snps_indels_svs.genotypes.tab | wc -l
echo "Create site removal list..."
#extract all those sites with same a0 allele but different a1 but only snps...so..not flipped alleles!..these need to be removed from our panel and from other panel if they are both singletons!
awk 'NR==FNR{b[\$2,\$4]=\$0;next;}{if(b[\$2,\$4]){split(b[\$2,\$4],x," ");if(\$5 != x[5] && length(\$5)==1) print b[\$2,\$4],\$4,\$5}}' \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.tab \${pop2_outpath}/\${pop2}.chr\${chr}.tab > \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching_ref.removelist

#now check if we need to remove those sites from the second panel because they are singletons in that panel
awk -v sAC=\${SAC} 'NR==FNR{b[\$2,\$8,\$9]=\$0;next;}{if(b[\$2,\$4,\$5]){if(\$7 == 1 || \$7==sAC) print \$0}}' \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching_ref.removelist \${pop2_outpath}/\${pop2}.chr\${chr}.tab > \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching_ref.\${pop2}_toremove

echo "Removing sites..."
#remove from vbi the mismatching snps
awk 'NR==FNR{b[\$2,\$4,\$5]=\$0;next;}!(b[\$2,\$4,\$5])' \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching_ref.removelist \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.tab > \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.purged.tab

#remove from the other panel the mismatching singleton snps
awk 'NR==FNR{b[\$2,\$4,\$5]=\$0;next;}!(b[\$2,\$4,\$5])' \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching_ref.\${pop2}_toremove \${pop2_outpath}/\${pop2}.chr\${chr}.tab > \${pop2_outpath}/\${pop2}.chr\${chr}.purged.tab

echo "Check singleton overlap..."
#extract singletons sites from sequences
awk '{if(\$7 == 1 || \$7 == 219) print \$0 }' \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.purged.tab > \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.purged.singleton.tab

#create keeplist (remove list) of singleton that overlap(not overlap) with 1Kgp
awk 'NR==FNR{b[\$2,\$4,\$5]=\$0;next;}!(b[\$2,\$4,\$5])' \${pop2_outpath}/\${pop2}.chr\${chr}.purged.tab \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.purged.singleton.tab > \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching.singleton.\${pop1}_toremove

echo "Removing singletons..."
#remove singleton not overlapping with 1kg (other panel) from vbi files
awk 'NR==FNR{b[\$2,\$4,\$5]=\$0;next;}!(b[\$2,\$4,\$5])' \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching.singleton.\${pop1}_toremove \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.purged.tab > \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.purged.sigremoved.tab

#extract singletons from 1kg and remove all that not overlapping with vbi
awk -v sAC=\${SAC} '{if((\$7 == 1 || \$7 == sAC)&&(length(\$4)==1 && length(\$5)==1)) print \$0 }' \${pop2_outpath}/\${pop2}.chr\${chr}.purged.tab > \${pop2_outpath}/\${pop2}.chr\${chr}.purged.singleton.tab

#create remove list
awk 'NR==FNR{b[\$2,\$4,\$5]=\$0;next;}!(b[\$2,\$4,\$5])' \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.purged.sigremoved.tab \${pop2_outpath}/\${pop2}.chr\${chr}.purged.singleton.tab > \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching.singleton.\${pop2}_toremove

#remove all non common singletons from 1kg(other) panel
awk 'NR==FNR{b[\$2,\$4,\$5]=\$0;next;}!(b[\$2,\$4,\$5])' \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching.singleton.\${pop2}_toremove \${pop2_outpath}/\${pop2}.chr\${chr}.purged.tab > \${pop2_outpath}/\${pop2}.chr\${chr}.purged.sigremoved.tab

echo "Indels purification for \${pop1}..."
#remove all indels overlapping with snps (by position) and all indels overlapping in the same file (by position) -> this are due to error from mpileup!!
fgrep -v -w -f <(cut -f 2 -d " " \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.purged.sigremoved.tab) \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.tab | cut -f 2 -d " " | uniq -d > \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.exclusion_list

#create the indel keeplist
fgrep -v -w -f \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.exclusion_list \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.tab > \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.purged.tab

#create list of allele mismatching sites
awk 'NR==FNR{b[\$2,\$4]=\$0;next;}{if(b[\$2,\$4]){split(b[\$2,\$4],x," ");if(\$5 != x[5] && length(\$5)!=1) print b[\$2,\$4],\$4,\$5}}' \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.purged.tab \${pop2_outpath}/\${pop2}.chr\${chr}.purged.sigremoved.tab > \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching_ref.indels.\${pop1}_toremove

#create list of allele mismatch sites that are singleton in the bigger panel
awk -v sAC=\${SAC} 'NR==FNR{b[\$2,\$8,\$9]=\$0;next;}{if(b[\$2,\$4,\$5]){if(\$7 == 1 || \$7 == sAC) print \$0}}' \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching_ref.indels.\${pop1}_toremove \${pop2_outpath}/\${pop2}.chr\${chr}.purged.sigremoved.tab > \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching_ref.indels.\${pop2}_toremove

#remove allele mismatching sites from VBI panel
awk 'NR==FNR{b[\$2,\$4,\$5]=\$0;next;}!(b[\$2,\$4,\$5])' \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching_ref.indels.\${pop1}_toremove \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.purged.tab > \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.purged.cleaned.tab

#remove allele mismatch from other panel if they are singletons (if the file is not empty)
if [ -s \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching_ref.indels.\${pop2}_toremove ]
then
awk 'NR==FNR{b[\$2,\$4,\$5]=\$0;next;}!(b[\$2,\$4,\$5])' \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.overlap_matching_ref.indels.\${pop2}_toremove \${pop2_outpath}/\${pop2}.chr\${chr}.purged.sigremoved.tab > \${pop2_outpath}/\${pop2}.chr\${chr}.purged.sigremoved.cleaned.tab
else
cp \${pop2_outpath}/\${pop2}.chr\${chr}.purged.sigremoved.tab \${pop2_outpath}/\${pop2}.chr\${chr}.purged.sigremoved.cleaned.tab
fi

#extract singleton indels(useless...)
awk '{if(\$7 == 1 || \$7 == 219) print \$0 }' \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.purged.cleaned.tab > \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.purged.cleaned.singleton.tab

#For indel sites we want to keep only those overlapping between panels and those not singletons
#create list of removable sites from vbi panel
awk 'NR==FNR{b[\$2,\$4,\$5]=\$0;next;}{if (!(b[\$2,\$4,\$5]) && (\$7==1 || \$7==219)) print \$0}' \${pop2_outpath}/\${pop2}.chr\${chr}.purged.sigremoved.cleaned.tab \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.purged.cleaned.tab > \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.indels.\${pop1}_toremove

#create the list of 1kg removable sites
awk -v sAC=\${SAC} 'NR==FNR{b[\$2,\$4,\$5]=\$0;next;}{if (!(b[\$2,\$4,\$5]) && (\$7==1 || \$7==sAC)) print \$0}' \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.purged.cleaned.tab \${pop2_outpath}/\${pop2}.chr\${chr}.purged.sigremoved.cleaned.tab > \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.indels.\${pop2}_toremove

echo "Removing INDELS...."
#remove sites from VBI
awk 'NR==FNR{b[\$2,\$4,\$5]=\$0;next;}!(b[\$2,\$4,\$5])' \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.indels.\${pop1}_toremove \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.purged.cleaned.tab > \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.purged.cleaned.final.tab

#remove sites from 1kgp
awk 'NR==FNR{b[\$2,\$4,\$5]=\$0;next;}!(b[\$2,\$4,\$5])' \${merged_outpath}/chr\${chr}.\${pop1}_\${pop2}.indels.\${pop2}_toremove \${pop2_outpath}/\${pop2}.chr\${chr}.purged.sigremoved.cleaned.tab > \${pop2_outpath}/\${pop2}.chr\${chr}.purged.sigremoved.cleaned.final.tab


#############################################
#now extract the selected sites from vcf file
#for 1kgp
function build_panel_extract_template(){
cat << EOF1
#!/usr/local/bin/bash

(zgrep "^#" \${1};awk 'NR==FNR{b[\\\$2,\\\$4,\\\$5]=\\\$0;next;}(b[\\\$2,\\\$4,\\\$5])' \${2} <(tabix \${1} \${5})) | bgzip -c > \${3}
tabix -f -p vcf \${3}
####now we can use the new vcf files to create legend files and haps files!!
#
#Legend
(echo "id position a0 a1";tabix \${3} \${5} | awk -v chrom=\${5} '{if (\\\$3 == ".") print "chr"chrom":"\\\$2,\\\$2,\\\$4,\\\$5;else print \\\$3,\\\$2,\\\$4,\\\$5}' | tr "\\t" " ") | gzip -c > \${4}.legend.gz
#haps
vcf-query \${3} -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GTR]\n' | tr "|" " " | tr "/" " " | tr "\\t" " "| cut -f 6- -d " " | gzip -c > \${4}.haps.gz

EOF1
}

echo "Creating first extraction script..."
build_panel_extract_template \${second_panel_file_path} \${pop2_outpath}/\${pop2}.chr\${chr}.purged.sigremoved.cleaned.final.tab \${pop2_outpath}/\${pop2}.chr\${chr}.genotypes.final.vcf.gz \${pop2_outpath}/HAPS/chr\${chr} \${chr} > \${pop2_outpath}/\${pop2}.chr\${chr}_extract_sites.sh

# bsub -J "extract_sites" -o"%J_extract.log" -e"%J_extract.err" -M3000000 -R"select[mem>3000] rusage[mem=3000]" -q normal -- bash \${pop2_outpath}/\${pop2}.chr\${chr}_extract_sites.sh

#for VBI (read the extract_sites.sh script for the command used)
#since we know that for snps there are not overlapping sites, we can extract the vcf desired sites by position!!
#these lines are meant to go in an extract_sites.sh script
function build_vbi_extract_template(){
cat << EOF2

#!/usr/local/bin/bash
(zgrep "^#" \${1};tabix \${1} \${9} | fgrep -w -f <(cut -f 1,2 -d " " \${2} |tr " " "\t")) | bgzip -c > \${3}

(zgrep "^#" \${4};awk 'NR==FNR{b[\\\$2,\\\$4,\\\$5]=\\\$0;next;}(b[\\\$2,\\\$4,\\\$5])' \${5} <(tabix \${4} \${9})) | bgzip -c > \${6}

tabix -f -p vcf \${3}
tabix -f -p vcf \${6}
#put together snps and indels for VBI
(zgrep "^#" \${3};(tabix \${3} \${9};tabix \${6} \${9}) | sort -g -k2,2)|bgzip -c > \${7}
tabix -f -p vcf \${7}
####now we can use the new vcf files to create legend files and haps files!!
#
#extract sites and create legend file and haps file
#VBI legend
(echo "id position a0 a1";tabix \${7} \${9} | awk -v chrom=\${9} '{if (\\\$3 == ".") print "chr"chrom":"\\\$2,\\\$2,\\\$4,\\\$5;else print \\\$3,\\\$2,\\\$4,\\\$5}' | tr "\\t" " ") | gzip -c > \${8}.legend.gz
#VBI haps
vcf-query \${7} -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GTR]\n' | tr "|" " " | tr "/" " " | tr "\\t" " "| cut -f 6- -d " " | gzip -c > \${8}.haps.gz

EOF2
}

echo "Creating second extraction script..."
build_vbi_extract_template \${first_panel_file_path} \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.purged.sigremoved.tab \${pop1_outpath}/\${pop1}.SNPS.chr\${chr}.final.vcf.gz \${indels_file_path} \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.purged.cleaned.final.tab \${pop1_outpath}/\${pop1}.INDELS.chr\${chr}.final.vcf.gz \${pop1_outpath}/\${pop1}.ALL.chr\${chr}.final.vcf.gz \${pop1_outpath}/HAPS/chr\${chr} \${chr} > \${pop1_outpath}/\${pop1}.chr\${chr}_extract_sites.sh

# bsub -J "extract_sites" -o"%J_extract.log" -e"%J_extract.err" -M3000000 -R"select[mem>3000] rusage[mem=3000]" -q normal -- bash \${pop1_outpath}/\${pop1}.chr\${chr}_extract_sites.sh

function build_merging_command(){
cat << EOF3
#!/usr/local/bin/bash
#command to create chunks and merge panels
#chese the start and end for chunk creation
pop1_start=\\\`zcat \${1} | head -2| tail -1| cut -f 2 -d " "\\\`
pop1_end=\\\`zcat \${1} | tail -1| cut -f 2 -d " "\\\`
pop2_start=\\\`zcat \${2} | head -2| tail -1| cut -f 2 -d " "\\\`
pop2_end=\\\`zcat \${2} | tail -1| cut -f 2 -d " "\\\`

if [ \\\${pop1_start} -ge \\\${pop2_start} ]
then
  chunk_start=\\\${pop2_start}
else
  chunk_start=\\\${pop1_start}
fi

if [ \\\${pop1_end} -ge \\\${pop2_end} ]
then
  chunk_end=\\\${pop1_end}
else
  chunk_end=\\\${pop2_end}
fi

impute_imputation_merging_s1.sh foo/bar \${6} \${3} \${4} \${5} \\\${chunk_start} \\\${chunk_end}

EOF3
}


echo "Creating merging command script..."
build_merging_command \${pop1_outpath}/HAPS/chr\${chr}.legend.gz \${pop2_outpath}/HAPS/chr\${chr}.legend.gz \${pop1_outpath}/HAPS \${pop2_outpath}/HAPS \${merged_outpath}/HAPS \${chr}> \${merged_outpath}/\${pop1}_\${pop2}.chr\${chr}_merge.sh

EOF
}

#build the panel input script
#we need as inputs:
#pop 1 (vbi) indel vcf filename
#pop 1 (vbi) snps vcf filename
#pop 2 vcf filename
#pop 1 output path
#pop 2 output path
#merged output
#
for chr in {1..22}
do
	pop1_snps_file=`ls ${first_panel_path}/*.chr${chr}.*.vcf.gz`
	pop1_indels_file=`ls ${indels_files_path}/*.chr${chr}.*.gz`
	#pop2_vcf_file=`ls ${second_panel_path}/*.chr${chr}.*.vcf.gz`
	pop2_vcf_file=`ls ${second_panel_path}/${chr}.*.20130411.vcf.gz`
	first_panel_outpath=${output_basepath}/${pop1}_${pop2}/${pop1}/CHR${chr}
	second_panel_outpath=${output_basepath}/${pop1}_${pop2}/${pop2}/CHR${chr}
	merged_output=${output_basepath}/${pop1}_${pop2}/MERGED/CHR${chr}

	mkdir -p ${first_panel_outpath}/HAPS
	mkdir -p ${second_panel_outpath}/HAPS
	mkdir -p ${merged_output}/HAPS

	build_panel_input_template ${pop1_snps_file} ${pop1_indels_file} ${pop2_vcf_file} ${first_panel_outpath} ${second_panel_outpath} ${merged_output} ${chr} ${pop1} ${pop2} > ${output_basepath}/${pop1}_${pop2}/merge_chr${chr}.sh

	
done

#test






# #create the folder
# mkdir -p HAPS/CHR7



# #1kgp (other panel)

# #now we can launch impute to create the merged panel!!
# #USAGE:\n impute_imputation_merging_s1.sh <genotype files path> <chr> <reference 1 files path> <reference 2 files path> <output file path> <chunk file ref>\n
# #- <genotype files path> : path for genotypes files (can be a fake path if you are only merging panels)
# #- <chr> : chromosome number
# #- <reference 1 files path> : path for first reference files
# #- <reference 2 files path> : path for second reference files
# #- <output files path> : output path
# #- <chunk file ref> : path to the reference file to use to create the chunks\n

# impute_imputation_merging_s1.sh foo/bar ${chr} ${ref_1_path}/HAPS/CHR${chr} ${ref_2_path}/HAPS/CHR${chr} ${ref_out_path}/CHR${chr}/HAPS ${ref_chunk_path}


# #create legend file for indels
# #(echo "id position a0 a1";awk '{if ($3 == ".") print "chr7:"$2,$2,$4,$5;else print $3,$2,$4,$5}' chr7.indels.keeplist | tr "\t" " ") > chr7.indels.legend

# #awk script for check of singleton and not INDEL overlap
# #1) extract only sites in common with 1000G with allele match
# # awk 'NR==FNR{a[$2,$3,$4]=$0;next;}(a[$2,$3,$4])' chr7_singletons.legend <( zcat /lustre/scratch109/sanger/mc14/GENOTIPI/REF_PANEL/1000G/ALL_1000G_phase1integrated_v3_impute/chr7.legend.gz)

# #2) extract only sites with same position and alleles mismatch
# # awk 'NR==FNR{a[$2,$3,$4]=$0;next;}(a[$2,$3,$4])' chr7_singletons.legend <( zcat /lustre/scratch109/sanger/mc14/GENOTIPI/REF_PANEL/1000G/ALL_1000G_phase1integrated_v3_impute/chr7.legend.gz)

# #3) extract sites with same position but flipped alleles
# # awk 'NR==FNR{a[$2,$3,$4]=$0;next;}(a[$2,$4,$3])' chr7_singletons.legend <( zcat /lustre/scratch109/sanger/mc14/GENOTIPI/REF_PANEL/1000G/ALL_1000G_phase1integrated_v3_impute/chr7.legend.gz)

# #extract singleton in legend file
# # cat chr7.legend.gz | fgrep -w -f <( cut -f 2 chr7_singletons.tab ) > chr7_singletons.legend

# #obtain singletons overlapping 1000gp
# # zcat ../1000G/ALL_1000G_phase1integrated_v3_impute/chr7.legend.gz | fgrep -w -f <(cut -f 2 -d " " chr7_singletons.legend) > chr7_singletons.overlap1kgp

# # vbi_singletons <- read.table('chr7_singletons.legend',sep=" ",header=F)
# # vbi_1kgp_singletons <- read.table('chr7_singletons.overlap1kgp',sep=" ",header=F)

# # colnames(vbi_1kgp_singletons) <- c("id","position","a0","a1","afr.aaf","amr.aaf","asn.aaf","eur.aaf","afr.maf","amr.maf","asn.maf","eur.maf")
# # colnames(vbi_singletons) <- c("id","position","a0","a1")

# # common_singletons_am <- vbi_singletons[which(vbi_singletons$position == vbi_1kgp_singletons$position & vbi_singletons$a0 == vbi_1kgp_singletons$a0 & vbi_singletons$a1 == vbi_1kgp_singletons$a1),]
