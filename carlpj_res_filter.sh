#!/usr/local/bin/bash

################################################
#script for CARLANTINO project
################################################
#

file=$1
filename=`basename ${file}`
chr1=${filename%%.*}
chr=${chr1#chr}
#Extract only significative results from analyses
# zcat ${file} | awk '$10<0.05' | gzip -c > ${filename}.sig.gz
# zcat ${file} | awk '{if ($10<0.05) print $3,$4,$5,$6,$2}'| awk '{if ($3=="-" || length($3) < length($4)) print $1,$2,$2+1,$3"/"$4,"+",$5;else if ($4=="-" || length($3) > length($4)) print $1,$2,$2+length($3)-1,$3"/"$4,"+",$5;else print $1,$2,$2,$3"/"$4,"+",$5;}' | gzip -c > ${filename}.sig.vep.gz

# #Annotate with VEP
# zcat ${filename}.sig.vep.gz | /nfs/team151/software/ensembl-tools-release-79/scripts/variant_effect_predictor/variant_effect_predictor.pl --format ensembl --o ${filename}.vep.annotated.tab --merged --all_refseq --gmaf --maf_esp --maf_1kg --quiet --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin Condel,/software/vertres/bin-external/VEP_plugins/config/Condel/config/,b --symbol --force_overwrite --html --cache --dir /data/blastdb/Ensembl/vep
# zcat ${filename}.sig.vep.gz | /nfs/team151/software/ensembl-tools-release-79/scripts/variant_effect_predictor/variant_effect_predictor.pl --format ensembl --o ${filename}.vep.annotated.tab --gmaf --maf_esp --quiet --regulatory --ccds --protein --uniprot --database --sift b --polyphen b --plugin Condel,/software/vertres/bin-external/VEP_plugins/config/Condel/config/,b --symbol --force_overwrite --html
#extract different consequences data
for conseq in missense UTR_variant stop_gained
do
	(grep ^# ${filename}.vep.annotated.tab;fgrep $conseq ${filename}.vep.annotated.tab) | gzip -c > ${filename}.$conseq.vep.annotated.tab.gz

	vep_in=${filename}.$conseq.vep.annotated.tab
	( zcat $vep_in.gz | head -100 | grep "^#"; zcat $vep_in.gz | grep -v "^#" | awk '{FS="\t";OFS="\t"}{gsub(/ /,"_",$0); print}' | awk '{FS="\t";OFS="\t"}{gsub(/:/,"\t",$2); print}'| sort -k2,2d -k3,3n ) | bgzip -c > $vep_in.compact.gz; 
	tabix -f -s 2 -b 3 -e 3 $vep_in.compact.gz

	awk 'FNR==NR { a[$1]=$0; next } $1 in a { print $0,a[$1] }' <(zcat ${filename}.sig.gz | tr " " "\t") <(tabix $vep_in.compact.gz ${chr} ) | gzip -c > ${filename}.$conseq.vep.annotated.tab.compact.merged.tab.gz
done

# file=/lustre/scratch113/teams/soranzo/users/jh21/imputed/carl/gemma/TG/chr19.gemma.gz


/nfs/team151/software/ensembl-tools-release-79/scripts/variant_effect_predictor/variant_effect_predictor.pl --format ensembl --o ${filename}.vep.annotated.tab --merged --all_refseq --gmaf --maf_esp --maf_1kg --quiet --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin Condel,/software/vertres/bin-external/VEP_plugins/config/Condel/config/,b --symbol --force_overwrite --html --cache --dir /data/blastdb/Ensembl/vep