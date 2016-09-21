#!/usr/local/bin/bash
# Angus Wallace, 2012
# wallace.angus@gmail.com
#   This program is free software: you can redistribute it and/or modify   #
#   it under the terms of the GNU General Public License as published by   #
#   the Free Software Foundation, either version 3 of the License, or      #
#   (at your option) any later version.                                    #
#   This program is distributed in the hope that it will be useful,        #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
#   GNU General Public License for more details.                           #
#   You should have received a copy of the GNU General Public License      #
#   along with this program.  If not, see <http: www.gnu.org="" licenses="">.  #
# version history
# v0.1 Angus Wallace Aug 2012
# v1.0 Angus Wallace Sep 2012

# a program to create Manhattan and QQ plots for huge GWAS data sets

### ############## ###
### check software ###
### ############## ###
if [ `which gnuplot | wc -l` -eq 0 ]; then echo "ERROR: gnuplot not found. Please install."; exit 1; fi
if [ `gnuplot -V | grep -o -E "[0-9]\.[0-9]" | awk '{ if($1 < 4.4) print 1}' | wc -m` -eq 1 ]; then #check gnuplot version number
   echo "ERROR: please install gnuplot version 4.4 or newer."; exit 1; fi
if [ `which awk | wc -l` -eq 0 ]; then echo "ERROR: awk not found. Please install."; exit 1; fi
if [ `which sort | wc -l` -eq 0 ]; then echo "ERROR: sort not found. Please install."; exit 1; fi
if [ `which sed | wc -l` -eq 0 ]; then echo "ERROR: sed not found. Please install."; exit 1; fi

### ################# ###
### process arguments ###
### ################# ###

if [ $# -eq 0 ]; then echo -e "\nCreate QQ and Manhattan plots for sanity-checking GWAS results. Written by Angus Wallace (wallace.angus@gmail.com), \n Use (switches are case-sensitive):\n  -in=filename    (specify a filename of input data)\n  -cwd            (work in the current directory)\n  -workdir=dir    (specify a directory in which to work)\n  -pvalueCol=num  (specify the column number of pvals in the datafile [default=11])\n  -SNPCol=n       (specify the column number of SNP names in the datafile [default=2])\n  -refdat=path    (specify the location of reference data)\n  -refChromoCol=n (specify the chromasome column in the reference datafile [default=1])\n  -refSNPCol=n    (specify the SNP column in the reference datafile [default=2])\n  -refBPCol=n     (specify the base-pair column in the reference datafile [default=4])\n  -nolabels       (specify that the datafile has no column labels)\n  -macImputed     (specify that the input data are produced by the minimac imputation program.  A separate reference data file is not needed, but can be supplied for faster execution)\n\nIf the input data filename has a full path, and neither -cwd or -workdir are used, then the path to the input file will be the working directory. If no paths are specified the current directory will be used.\n\nExamples:\nRun fastplot on the output of Merlin for a particular trait, using the standard reference data set:\n  fastplot -in=trait_data.out -refdat=/scratch/GWAS/GeneralRelease/Release6/GWAS.bim\nRun fastplot on the output of minimac imputed 1000Genome data:\n  fastplot -in=subcort.dat -macImputed\nRun fastplot on some data with columns arranged differently to the defaults:\n  fastplot -in=trait_data.out -refdat=/scratch/GWAS/GeneralRelease/Release6/GWAS.bim -pvalueCol=3 -SNPCol=2 -nolabels refChromoCol=3\n\nThis program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program.  If not, see http://www.gnu.org/licenses.\n"; exit 1; fi
 
col_lab=1 #assume the data file has labels unless told otherwise
imputed=0 #assume the data are not imputed unless told otherwise
subsample=0
sort_mem=50M
for arg in "$@" ; do
   case $arg in
      -in=*) # input file name
         in_dat=`echo $arg | sed 's/[-a-zA-Z0-9]*=//'`
         ;;
      -cwd) #work in the current directory
         workdir=`pwd`
         ;;
      -workdir=*) #specify a directory in which to work
         workdir=`echo $arg | sed 's/[-a-zA-Z0-9]*=//'`
         ;;
      -pvalueCol=*) #specify which column has the pvalues
         pvalueCol=`echo $arg | sed 's/[-a-zA-Z0-9]*=//'`
         ;;
      -SNPCol=*) #specify which column has the SNP names
         SNPCol=`echo $arg | sed 's/[-a-zA-Z0-9]*=//'`
         ;;
      -refdat=*) #specify reference data
         ref_dat=`echo $arg | sed 's/[-a-zA-Z0-9]*=//'`
         ;;
      -refChromoCol=*) #specify reference data chromasome col
         refChromoCol=`echo $arg | sed 's/[-a-zA-Z0-9]*=//'`
         ;;
      -refSNPCol=*) #specify reference data SNP col
         refSNPCol=`echo $arg | sed 's/[-a-zA-Z0-9]*=//'`
         ;;
      -refBPCol=*) #specify reference data's base-pair col
         refBPCol=`echo $arg | sed 's/[-a-zA-Z0-9]*=//'`
         ;;
      -nolabels) #line one has column labels
         col_lab=0
         ;;
      -macImputed) #line one has column labels
         imputed=1
         ;;
      -sortmem=*) # specify how much RAM is allocated to sort commands
         sort_mem=`echo $arg | sed 's/[-a-zA-Z0-9]*=//'`
         ;;
#      -subsample) #line one has column labels
#         subsample=1
#         ;;
   esac
done 
 
# if there's no path in the filename, then use the cwd
if [ `echo ${workdir} | wc -m` -eq 1 ]; then
   workdir=`echo ${in_dat} | grep -E -o ".*/"`
   if [ `echo ${workdir} | wc -m` -eq 1 ]; then
      workdir=`pwd`; fi; fi
cd $workdir
 
if [ `echo ${in_dat} | wc -m` -eq 1 ]; then
   echo "ERROR: no input file specified."
   exit 1
elif [ `ls ${workdir}/${in_dat} | wc -m` -eq 0 ]; then
   echo "ERROR: input file not found."
   exit 1
fi
  
 
 
### ######## ###
### defaults ###
### ######## ###
 
if [ `echo ${pvalueCol} | wc -m` -eq 1 ]; then pvalueCol=11; fi
if [ `echo ${SNPCol} | wc -m` -eq 1 ]; then SNPCol=2; fi
if [ `echo ${refChromoCol} | wc -m` -eq 1 ]; then refChromoCol=1; fi
if [ `echo ${refSNPCol} | wc -m` -eq 1 ]; then refSNPCol=2; fi
if [ `echo ${refBPCol} | wc -m` -eq 1 ]; then refBPCol=4; fi
 
if [ `echo ${ref_dat} | wc -m` -eq 1 ]; then
   if [ $imputed -eq 0 ]; then
      echo "ERROR: no reference file specified. try -refdat=/scratch/GWAS/GeneralRelease/Release6/GWAS.bim"
      exit 1; 
   else #then the data were inputed with minimac
      # create reference data from input data SNP names
      awk " \$${pvalueCol} != \"-\" {sub(\".*:\",\"\",\$2); print \$1, \$1\":\"\$2, \$2}" $in_dat | sort  -S ${sort_mem} -g > fastplot_refdata.dat
# possible way of specifying different columns for imputed data.
#      awk " \$${pvalueCol} != \"-\" {\$0= \$${SNPCol}\" \"\$${SNPCol}\" \"\$${SNPCol}; sub(\":.*\",\"\",\$1);sub(\".*:\",\"\",\$3); print \$0}" $in_dat | sort  -S ${sort_mem} -g > fastplot_refdata.dat
      ref_dat=${workdir}/fastplot_refdata.dat
      refBPCol=3
   fi
fi
 
### ####################################### ###
### create the Manhattan data and plot them ###
### ####################################### ###
 
# from ref_dat, create a list of cumulative SNP counts across chromasomes
rm -f fastplot_chr_len.dat
for i in `cat ${ref_dat} | awk '{print $1}' | sort -S ${sort_mem} | uniq | sort -g | uniq`; do
   awk "BEGIN {max = 0} {if (\$1==${i} && \$${refBPCol}>max) max=\$${refBPCol}} END {print max}" ${ref_dat} >> fastplot_chr_len.dat
done
 
echo "1 0" > abs_BP.dat
awk 'BEGIN{sum=0}{sum += $1; print (NR+1), sum}' fastplot_chr_len.dat | head -n 21 >> abs_BP.dat
 
chr_len=abs_BP.dat
 
if [ $imputed -eq 0 ]; then #are we dealing with the imputed output of mac?
   # prune input data to pvals column, remove first line, sort
   if [ $col_lab -eq 1 ]; then
      awk "NR > 1 && \$${pvalueCol} != \"-\" {print \$${SNPCol}, \$${pvalueCol} } " ${in_dat} | sort  -S ${sort_mem} > ${workdir}/manh_in_sorted.dat
   else
      awk "\$${pvalueCol} != \"-\" {print \$${SNPCol}, \$${pvalueCol} } " ${in_dat} | sort  -S ${sort_mem} > ${workdir}/manh_in_sorted.dat
   fi
   # get the base-pairs of the SNPs, translate to absolute BPs (not relative to chr)
   awk "{print \$${refChromoCol}, \$${refSNPCol}, \$${refBPCol} }" ${ref_dat} > ${workdir}/ref_dat.dat
   join ${workdir}/ref_dat.dat $chr_len | awk {'print $2, ($3 + $4), ($1 % 2)+1'} | sort  -S ${sort_mem} > ${workdir}/step_1_all.dat
 
   #merge the two files crop to pval/BP and take the -log10
   join ${workdir}/manh_in_sorted.dat ${workdir}/step_1_all.dat | awk {'print $3, -log($2) / log(10), $4'} > ${workdir}/all_dat.dat
   ### clean up ###
   rm ${workdir}/step_1_all.dat ${workdir}/manh_in_sorted.dat 
else
   # then we're using minimac (13M) data. Do as much as possible in one line for speed.
   awk "NR==FNR{a[\$1]=\$2;next;}{  {sub(\$1,a[\$1]\" \"\$1) sub(\".*:\",\"\",\$3)}; print \$1+\$3, (-log(\$12 )/2.302585), (\$2 % 2)+1 }" abs_BP.dat ${in_dat} > all_dat.dat
#   fi
fi
 
# create the plot
if [ `gnuplot -V | grep -o -E "[0-9]\.[0-9]" | awk '{ if($1 < 4.4) print 1}' | wc -m` -ge 1 ]; then
   echo -e "\nWARNING: gnuplot 4.4 or newer is needed for prettier Manhattan plot output, colour-coded by chromosome\n"
gnuplot -persist <<plot set="" term="" png="" size="" 1280,600="" unset="" colorbox="" output="" "manhattan_output_simple.png"="" ylabel="" "-log10(pvalue)"="" xlabel="" "absolute="" bp"="" title="" "manhattan="" plot"="" plot="" "all_dat.dat"="" u="" 1:2="" with="" points="" else="" gnuplot="" -persist="" <<plot="" "manhattan_output.png"="" palette="" model="" rgb="" defined="" (="" 1="" 'red',="" 2="" 'green'="" )="" 1:2:3="" notitle="" pt="" 6="" fi="" ###="" clean="" up="" rm="" all_dat.dat="" ################################="" create="" the="" qq="" data="" and="" them="" #="" prune="" to="" pvals="" column,="" remove="" first="" line,="" rows="" without="" pvalue="" estimate,="" sort="" using="" scientific="" notation,="" take="" log10="" awk="" "{if="" (\$${pvaluecol}="" !="\"-\")" {print="" (-log(\$${pvaluecol})="" 2.302585)="" }}"="" ${in_dat}="" |="" -s="" ${sort_mem}="" -g=""> ${workdir}/qq_in_sorted.dat
 
# how many data are there?
dat_len=`cat ${workdir}/qq_in_sorted.dat | wc -l`
 
# generate chi_squared data of the same length
seq ${dat_len} -1 1 | awk {'print -log($1/'$dat_len')/ 2.302585'}  > ${workdir}/chi2_temp.dat # ${workdir}/ints.dat
 
# combine into the data to be plotted
paste ${workdir}/chi2_temp.dat ${workdir}/qq_in_sorted.dat > ${workdir}/chi2_plot_data.dat
paste ${workdir}/chi2_temp.dat ${workdir}/chi2_temp.dat > ${workdir}/x_eq_y.dat
 
# create the QQ plot
gnuplot -persist <<plot set="" term="" png="" output="" "qq_output.png"="" ylabel="" "observed"="" xlabel="" "expected"="" title="" "qq-plot"="" plot="" "chi2_plot_data.dat"="" with="" points="" notitle,="" "x_eq_y.dat"="" lines="" notitle="" key="" left="" top="" ###="" calculate="" lambda="" #="" trim="" both="" the="" observed="" and="" expected="" data="" at="" median,="" find="" mean="" of="" remaining="" data,="" ratio="" between="" is="" trim_len="$(($dat_len" 2="" ))="" awk="" "begin{sum1="0;" sum2="0}" nr<\$trim_len="" {sum1="" +="\$1" ;sum2="" printf\"%0.10f="" %0.10f\\n\",="" sum1,="" sum2}"="" ${workdir}="" chi2_plot_data.dat="" |="" tail="" -n="" 1="" '{printf"%0.10f\n",="" ($2="" $1)="" }'=""> ${workdir}/lambda_val.txt
 
echo "QQ lambda estimate: " `cat ${workdir}/lambda_val.txt`
 
 
### clean up ###
rm chi2_plot_data.dat x_eq_y.dat chi2_temp.dat qq_in_sorted.dat
 
