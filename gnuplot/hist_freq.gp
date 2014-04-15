#to plot histogram/density like 
clear
reset
set key off
set border 3

# Add a vertical dotted line at x=0 
set yzeroaxis
set boxwidth 0.0005 absolute
set style fill solid 1.0 noborder
bin_width = 0.001;

bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
plot "./chr20.frq" using (rounded($5)):(5) smooth frequency with boxes
