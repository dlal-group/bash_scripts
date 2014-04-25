#!/usr/local/bin/bash
#to create a gnuplot script to plot histogram/density like 
#Args:
#$1=filename
#$2=columns

function create_gnuplot_script() {
cat << EOF
clear
reset
filename="`echo "$1"`"
set terminal png
set output "`echo $1.png`"
set key off
# set border 3

# Add a vertical dotted line at x=0 
# set yzeroaxis
set boxwidth 0.0005 absolute
# set style fill solid 1.0 noborder
bin_width = 0.001;

bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
plot filename using (rounded(\$$2)):($2) smooth frequency with boxes

EOF
}

create_gnuplot_script $1 $2 > $1.gnuplot

gnuplot $1.gnuplot