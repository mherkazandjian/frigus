unset label
set log
set terminal pdf monochrome solid font 'Helvetica,14' size 16cm,12cm
set output "h2_cf_updated.pdf"
set format y "10^{%T}" 
set xlabel "T_g"
set ylabel "cooling function"
set mytics 10
set key bottom right

plot 'cooling_nc=1.0E+06' u 2:3 w l ls 1 ti 'nc = 1.0E+06', \
'cooling_nc=1.0E+14' u 2:4 w p ti 'fit by Glover and Abel 2008'
