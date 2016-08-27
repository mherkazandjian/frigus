unset label
set log
set output "cooling_function.eps"
set term postscript eps monochrome enhanced "Helvetica" 20
set format y "10^{%T}" 
set xlabel "T_g"
set ylabel "cooling function"
#set xrange [2.e3:2e1]
#set yrange [1e-20:1.e-5]
set key bottom right
plot 'cooling.txt' u 2:3  ti 'present calculation', 'cooling.txt' u 2:4  ti 'glover'

