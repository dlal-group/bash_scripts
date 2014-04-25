#!/usr/local/bin/bash
#to create a gnuplot script to plot histogram/density like 
#Args:
#$1=filename
#$2=expected value
#$3=observed value


function create_gnuplot_script() {
cat << EOF
#Gnuplot script for qqplot and manhattan plot from data in a file
clear
reset
filename="`echo "$1"`"
set terminal png
set output "`basename $1\_qqplot.png`"
set key off
set title "QQ plot"
unset autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
expmax = `tail -n+2 $1 | awk '{if(min==""){min=max=\$7}; if(\$7>max) {max=\$7}; if(\$7< min) {min=\$7}; total+=\$7; count+=1}END{print max}'`
obsmax = `tail -n+2 $1 | awk '{if(min==""){min=max=\$5}; if(\$5>max) {max=\$5}; if(\$5< min) {min=\$5}; total+=\$5; count+=1}END{print max}'`
set xtic 0,1,expmax                     # set xtics automatically
set ytic 0,10,obsmax                         # set ytics automatically
set xlabel "Expected -log[10](p)"
set ylabel "Observed -log[10](p)"
set xr [0.0:expmax]
set yr [0:obsmax]
#to plot the expected line we need to replicate this:
#expmax <- trunc(max(logexppval))+1

plot filename using $2:$3 with points, \
	filename using $2:$2 with lines

EOF
}

filename=`basename $1`
create_gnuplot_script $1 $2 $3 > ${filename}.gnuplot

#gnuplot $1.gnuplot