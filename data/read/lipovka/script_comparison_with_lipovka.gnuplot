set log
set output "comparison_with_lipovka_tab.eps"
set term postscript eps color enhanced "Helvetica" 20
set format y "10^{%T}" 
set xlabel "T_g"
set ylabel "cooling function"
#set xrange [2.e3:2e1]
#set yrange [1e-20:1.e-5]
set key bottom right
plot 'cooling_function_lipovka.dat' u (10**$1):(1e-7*10**$3) ti 'lipovka data n_H = 10^7 m-3', 'cooling_nc=1.0E+07' u 2:3 ti '', 'cooling_function_lipovka.dat' u (10**$1):(1e-7*10**$4) ti 'lipovka data n_H = 10^8 m-3', 'cooling_nc=1.0E+08' u 2:3 ti '', 'cooling_function_lipovka.dat' u (10**$1):(1e-7*10**$5) ti 'lipovka data n_H = 10^9 m-3', 'cooling_nc=1.0E+09' u 2:3 ti '', 'cooling_function_lipovka.dat' u (10**$1):(1e-7*10**$6) ti 'lipovka data n_H = 10^{10} m-3', 'cooling_nc=1.0E+10' u 2:3 ti '', 'cooling_function_lipovka.dat' u (10**$1):(1e-7*10**$7) ti 'lipovka data n_H = 10^{11} m-3', 'cooling_nc=1.0E+11' u 2:3 ti '', 'cooling_function_lipovka.dat' u (10**$1):(1e-7*10**$9) ti 'lipovka data n_H = 10^{14} m-3', 'cooling_nc=1.0E+14' u 2:3 ti ''

