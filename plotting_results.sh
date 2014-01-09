#!/usr/local/bin/bash
#ARGS:
#$1=trait
#$2=filepath

trait=$1
res_path=$2
pos_con_list=$3

#IFS=$'\n'
#set $(cat parameter_file.txt)

# RUN GenABEL.R
#bsub -J "G_$2_$4" -o GenABEL.log  -M7000000 -R"select[mem>7000] rusage[mem=7000]" \
#-q long R CMD BATCH '--args '$6' '$4' '$8' '${10}' '${12}' '$2'' /nfs/team151/GWAS/IMPUTATION-1KG/VALBORBERA/GWAS/GenABEL.R

# RUN ProbABEL
#for i in $(seq 1 22)
#do
#bsub -J "P_$2_$4_chr${i}" -w "ended(G_$2_$4)" -o chr${i}.log -M7000000 -R"select[mem>7000] \
#rusage[mem=7000]" -q basement /nfs/users/nfs_g/gp5/ProbABEL/bin/palinear --mmscore varcovar.mat --pheno res.pheno \
#--chrom ${i} --map ${14}/${i}/mach.${i}.machlegend \
#--info ${16}/${i}/mach.${i}.machinfo \
#--dose ${18}/output${i}.gen.dose.fvi --out chr${i}.palinear
#done

# RUN merge.R
#bsub -o merge.log -J "merge_$2_$4" -w "ended(P_$2_$4_chr*)" \
#-M10000000 -R"select[mem>10000] rusage[mem=10000]" -q hugemem \
#R CMD BATCH '--args '$2' '$4'' /nfs/team151/GWAS/IMPUTATION-1KG/VALBORBERA/GWAS/merge.R

# RUN exclude.pl
#bsub -o exclude.log -J "PERL_$2_$4" -w "ended(merge_$2_$4)" \
#-M6000000 -R"select[mem>6000] rusage[mem=6000]" -q long \
#perl /nfs/users/nfs_a/ar10/exclude.pl all.palinear.out all.palinear.filtered
# function build_gnuplot_script(){
# cat << EOF
# #Gnuplot script for qqplot and manhattan plot from data in a file
# set terminal png
# set output "$1_qqplot.png"
# #set autoscale                        # scale axes automatically
# unset autoscale                        # scale axes automatically
# unset log                              # remove any log-scaling
# unset label                            # remove any previous labels
# lexpmax = `fgrep -v inf $1 | awk -v col=$2 '{if(min==""){min=max=\$(col)}; if(\$(col)>max) {max=$(col)}; if($(col)< min) {min=$(col)}; total+=$(col); count+=1}END{print max}'`
# lobsmax = `fgrep -v inf $1 | awk -v col=$3 '{if(min==""){min=max=\$(col)}; if(\$(col)>max) {max=$(col)}; if($(col)< min) {min=$(col)}; total+=$(col); count+=1}END{print max}'`
# lobs=$3
# lexp=$2
# set xtic 0,1,lexpmax                     # set xtics automatically
# set ytic 0,10,lobsmax                         # set ytics automatically
# set key off
# set title "QQ plot $1"
# set xlabel "Expected -log[10](p)"
# set ylabel "Observed -log[10](p)"
# set xr [0.0:lexpmax]
# set yr [0:lobsmax]
# #to plot the expected line we need to replicate this:
# #expmax <- trunc(max(logexppval))+1

# plot "$1" using lexp:lobs with points, \
#         "$1" using lexp:lexp with lines
# EOF
# }

# file=$1
# lexp=$2
# lobs=$3

# # RUN FORMATTING AND PLOTTING STEPS
# build_gnuplot_script ${file} ${lexp} ${lobs} > ${file}_gnuplot.gp

# gnuplot ${file}_gnuplot.gp

#for chr in {1..22} X
#do
# res_file=${res_path}/${trait}.all.tab.assoc.1e-1.to_plot
#we need to know how many sites in total we have before the downsampling...
# tot_sites=`wc -l ${res_path}/${trait}.*.result.maf_info_hwe_filtered.join | grep "total"| awk '{print $1}'`
res_file=${res_path}/${trait}.all.result.maf_info_hwe_filtered.join.plot

# R CMD BATCH '--args '${trait}' '${res_file}' '${pos_con_list}' '${tot_sites}'' /nfs/users/nfs_m/mc14/Work/r_scripts/replica_plotter.r
R CMD BATCH '--args '${trait}' '${res_file}' '${pos_con_list}'' /nfs/users/nfs_m/mc14/Work/r_scripts/replica_plotter.r ${res_path}/${trait}.Rout

#done

