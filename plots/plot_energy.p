set terminal wxt size 900,700 enhanced font 'Verdana,10' persist
set multiplot layout 2,1 rowsfirst


#gráfico a
set label 1  
set xlabel 'theta'

set ylabel 'dE mu+'

file(n) = sprintf("energy_mu_plus_data/%d0.000000.dat",n)

plot for [i=1:10] file(i) u 2:1 title sprintf("%d0 GeV",i)




#grafico 2

set label 2 
set xlabel 'theta'
set ylabel 'dE mu-' 
file(n) = sprintf("energy_mu_minus_data/%d0.000000.dat",n)

plot for [i=1:10] file(i) u 2:1 title sprintf("%d0 GeV",i)
