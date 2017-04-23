unset label
set log
set output "cf_lipovka.eps"
set term postscript eps monochrome enhanced "Helvetica" 20
set format y "10^{%T}" 
set xlabel "T_g"
set ylabel "cooling function"
set mytics 10
#set xrange [2.e1:2e3]
#set yrange [1e-20:1.e-5]
set key bottom right
#plot 'cooling_nc=1.0E+06' u 2:3 w l ti 'nc = 1.0E+06', \
#'cooling_nc=1.0E+07' u 2:3 w l ti 'nc = 1.0E+07',      \
#'cooling_nc=1.0E+08' u 2:3 w l ti 'nc = 1.0E+08',      \
#'cooling_nc=1.0E+09' u 2:3 w l ti 'nc = 1.0E+09',      \
#'cooling_nc=1.0E+10' u 2:3 w l ti 'nc = 1.0E+10',      \
#'cooling_nc=1.0E+11' u 2:3 w l ti 'nc = 1.0E+11',      \
#'cooling_nc=1.0E+12' u 2:3 w l ti 'nc = 1.0E+12',      \
#'cooling_nc=1.0E+13' u 2:3 w l ti 'nc = 1.0E+13',      \
#'cooling_nc=1.0E+14' u 2:3 w l ti 'nc = 1.0E+14',      \
#'cooling_nc=1.0E+14' u 2:4 w l ti 'fit by Glover and Abel 2008',      \
#'mher_cooling.txt' u 1:2 w p ti 'result by Mher @ T_r = 3000 K', \

plot 'cooling_nc=1.0E+14' u 2:3 w l ti 'nc = 1.0E+14',\
'cooling_nc=1.0E+14' u 2:4 w l ti 'fit by Lipovka',\
'cooling_nc=1.0E+06' u 2:3 w l ti 'nc = 1.0E+06',\
'cooling_nc=1.0E+06' u 2:4 w l ti 'fit by Lipovka'
