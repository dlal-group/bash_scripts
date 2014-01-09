#!/usr/local/bin/bash
#Launch qc r script via LSF
#TODO: use getopt to read specified options...
function create_plotter_file() {
cat << EOF
#script to plot all chr files together
rm(list=ls())
require(gplots)

#Set arguments to use from command line
# 
all_VBI <- NULL
all_TSI <- NULL
all_vbi_snps <- NULL
prefix_vbi <- "$1"
prefix_alternative <- "$2"
suffix <- "$3"
prefix_vbi_snps <- "$4"
suffix_vbi_snps <- "$5"

for(i in 1:22){
#	vbi_file_path <- paste(prefix_vbi,i,suffix,sep="")
#	alternative_file_path <- paste(prefix_alternative,i,suffix,sep="")
	#READ ALL DATA
#	current_VBI_snps <- read.table(vbi_file_path, sep='\t',header=T)
#	current_alternative_snps <- read.table(alternative_file_path, sep='\t',header=T)
#	all_VBI <- rbind(all_VBI,current_VBI_snps)
#	all_EUR <- rbind(all_EUR,current_alternative_snps)

	vbi_snps_file_path <- paste(prefix_vbi_snps,i,suffix_vbi_snps,sep="")
	all_current_vbi_snps <- read.table(vbi_snps_file_path, sep='\t',header=T)
	all_vbi_snps <- rbind(all_vbi_snps,all_current_vbi_snps)
	rm(all_current_vbi_snps)
	gc()
}

#calculate the maf
vbi_af_over <- all_vbi_snps[which(all_vbi_snps\$AF > 0.5),]
all_vbi_snps <- all_vbi_snps[-which(all_vbi_snps\$AF > 0.5),]
#now transform the column value in vbi_af_over
vbi_af_over\$AF <- (1-vbi_af_over\$AF)
#now rebuilt the complete dataset again
all_vbi_snps <- rbind(all_vbi_snps,vbi_af_over)

#plot
#jpeg("maf_density_overlap_ALL.jpg", width=1000, height=1000)
#	plot(density(all_VBI\$VBI_AF), main="Maf distribution",xlab="MAF",ylab="Maf density distribution",col=c("red"))
#	lines(density(all_EUR\$EUR_AF),col=c("blue"))
#	smartlegend(x="right",y="top", inset = 0,c("VBI", "EUR"),fill = c("red", "blue"))
#dev.off()
#jpeg("maf_barplot_overlap_ALL_EUR.jpg", width=1000, height=1000)
#	plot(as.factor(all_EUR\$EUR_AF), main="Maf by classes in EUR",xlab="MAF",ylab="Snps number",col=c("Blue"))
#dev.off()
#jpeg("maf_barplot_overlap_ALL_VBI.jpg", width=1000, height=1000)
#	plot(as.factor(all_VBI\$VBI_AF), main="Maf by classes in VBI",xlab="MAF",ylab="Snps number",col=c("Red"))
#dev.off()
jpeg("maf_density_ALL_SNPS_VBI.jpg", width=1000, height=1000)
	plot(density(as.numeric(as.character(all_vbi_snps\$AF))), main="Maf distribution",xlab="MAF",ylab="Maf density distribution",col=c("red"))
dev.off()
jpeg("maf_barplot_ALL_SNPS_VBI.jpg", width=1000, height=1000)
	plot(as.factor(all_vbi_snps\$AF), main="Maf by classes in VBI",xlab="MAF",ylab="Snps number",col=c("Red"))
dev.off()
q(save="no")
EOF
}

#FIXME:make this smart:we need to create the directory tree only if there are / in the string
mkdir -p QC_OUT
cd QC_OUT/
mkdir -p SNPS_FILES
mkdir -p SNPS_FILES/MAF_FILES
mkdir -p SNPS_FILES/VBI_SNPS_NO_MULTI
mkdir -p SNPS_FILES/VQSLOD
mkdir -p PDF
mkdir -p JPG
mkdir -p MONO
mkdir -p MULTIALLELIC
mkdir -p MULTIALLELIC/SUMMARY
mkdir -p OVERLAP
mkdir -p OVERLAP/COMPLETE_OVERLAP
mkdir -p OVERLAP/DIFFS
mkdir -p OVERLAP/PLOT_TABLES
mkdir -p SUMMARIES

echo "Starting QC script..."
#-J job_name_spec[index |  start_index-end_index:step,]
vbi_file_path='/nfs/users/nfs_m/mc14/lustre_home/GENOTIPI/COMPARISON/VBSEQ_QC/VBI/VBI.23chr.snps.reann.vcf'
population_file_path='/nfs/users/nfs_m/mc14/lustre_home/GENOTIPI/COMPARISON/VBSEQ_QC/TSI/TSI.23chr.snps.no_AC_zero.vcf'

#extract header from file
vbi_head=`head -1 $vbi_file_path`
population_head=`head -1 $population_file_path`

#split file by chr and submit jobs by chr
for chr in {1..22}
do
	echo "Extract chr${chr} info..."
	cat <(echo $vbi_head) <(egrep "^${chr}	" $vbi_file_path) > 'VBI_chr'${chr}'.snps.csv'
	sed -i 's/ /\t/g' 'VBI_chr'${chr}'.snps.csv'
	cat <(echo $population_head) <(egrep "^${chr}	" $population_file_path) > 'TSI_chr'${chr}'.snps.csv'
	sed -i 's/ /\t/g' 'TSI_chr'${chr}'.snps.csv'
	mv *.snps.csv SNPS_FILES/
	bsub -J "QC_script_chr_${chr}" -o "%J_chr${chr}_out.log"  -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
	-q basement R CMD BATCH "--args SNPS_FILES/VBI_chr${chr}.snps.csv SNPS_FILES/TSI_chr${chr}.snps.csv ${chr}" /nfs/users/nfs_m/mc14/lustre_home/GENOTIPI/COMPARISON/VBSEQ_QC/scripts/esgi_vbseq_qc_TSI.r 
done

#bsub -J "QC_script_TSI_chr_[1-22]" -o "%J_chr%I_out.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
#-q basement 'R CMD BATCH "--args SNPS_FILES/VBI_chr"${LSB_JOBINDEX}".snps.csv SNPS_FILES/TSI_chr"${LSB_JOBINDEX}".snps.csv "${LSB_JOBINDEX}"" /nfs/users/nfs_m/mc14/lustre_home/GENOTIPI/COMPARISON/VBSEQ_QC/scripts/esgi_vbseq_qc_TSI.r'

#create plotter file runtime
#set prefixes and suffixes for files to read..
#in the script we want sthing like this...
#	vbi_file_path <- paste('VBI_chr',i,'_plot_table_MAF.csv',sep="")
#	alternative_file_path <- paste('EUR_chr',i,'_plot_table_MAF.csv',sep="")

create_plotter_file "VBI_chr" "TSI_chr" "_plot_table_MAF.csv" "VBI_snps_chr" "_table.vcs"> /nfs/users/nfs_m/mc14/lustre_home/GENOTIPI/COMPARISON/VBSEQ_QC/scripts/plotter_TSI.R

#extract jobarray id
#jarrayid=`bjobs -q basement -A | cut -f 1 -d " " | tail -n 1`

#bsub -J "QC_plot_TSI" -w "numdone($jarrayid,22)" -o "QC_plot_out.log" -M16000000 -R"select[mem>16000] rusage[mem=16000]" \
#-q basement R CMD BATCH /nfs/users/nfs_m/mc14/lustre_home/GENOTIPI/COMPARISON/VBSEQ_QC/scripts/plotter_TSI.R

#bsub -J "move_job_TSI" -w "ended(QC_plot_TSI)" -o "move_job_out.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
#-q basement "mv *.pdf PDF/;mv *.jpg JPG/;mv multiallelic_chr*_table.csv MULTIALLELIC/;mv *mono_RR_chr* MONO/;mv *mono_AA_chr* MONO/;mv *plot_table_MAF* OVERLAP/PLOT_TABLES/;mv VBI_TSI_diff_rsIDs* OVERLAP/DIFFS/;mv VBI_TSI_all_overlap* OVERLAP/COMPLETE_OVERLAP/;mv *_summary* SUMMARIES/;mv VBI_snps_chr*_table.vcs SNPS_FILES/VBI_SNPS_NO_MULTI/"
