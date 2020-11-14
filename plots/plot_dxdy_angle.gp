#set terminal wxt size 1000,900 enhanced font 'Verdana,10' persist

set terminal postscript enhanced color solid  "Times" 14
set output 'dxdy.eps'

set multiplot layout 2,1 rowsfirst


#gráfico 1
set label 1
set key top left
set xlabel 'Energy(GeV)'
set ylabel '{/Symbol D}x'
set xrange[0:11]
set yrange[-0.15:1]
plot "data_mu_plus_position_1.4.dat" u 1:7 title "1.4", "data_mu_plus_position_1.9.dat" u 1:7 title "1.9", "data_mu_plus_position_2.2.dat" u 1:7 title "2.2", "data_mu_plus_position_2.7.dat" u 1:7 title "2.7", "data_mu_plus_position_3.0.dat" u 1:7 title "3.0", "data_mu_plus_position_3.5.dat" u 1:7 title "3.5", "data_mu_plus_position_4.3.dat" u 1:7 title "4.3", "data_mu_plus_position_5.1.dat" u 1:7 title "5.1", "data_mu_plus_position_5.7.dat" u 1:7 title "5.7", "data_mu_plus_position_6.3.dat" u 1:7 title "6.3", "data_mu_plus_position_7.1.dat" u 1:7 title "7.1", "data_mu_plus_position_8.1.dat" u 1:7 title "8.1", "data_mu_plus_position_8.7.dat" u 1:7 title "8.7", "data_mu_plus_position_9.0.dat" u 1:7 title "9.0", "data_mu_plus_position_9.4.dat" u 1:7 title "9.4", "data_mu_plus_position_9.7.dat" u 1:7 title "9.7", "data_mu_plus_position_9.9.dat" u 1:7 title "9.9"

#grafico 2
set key top right
set label 1 
set ylabel '{/Symbol D}y'
set xrange[0:11]
set yrange[-0.4:1]
plot "data_mu_plus_position_1.4.dat" u 1:8 title "1.4", "data_mu_plus_position_1.9.dat" u 1:8 title "1.9", "data_mu_plus_position_2.2.dat" u 1:8 title "2.2", "data_mu_plus_position_2.7.dat" u 1:8 title "2.7", "data_mu_plus_position_3.0.dat" u 1:8 title "3.0", "data_mu_plus_position_3.5.dat" u 1:8 title "3.5", "data_mu_plus_position_4.3.dat" u 1:8 title "4.3", "data_mu_plus_position_5.1.dat" u 1:8 title "5.1", "data_mu_plus_position_5.7.dat" u 1:8 title "5.7", "data_mu_plus_position_6.3.dat" u 1:8 title "6.3", "data_mu_plus_position_7.1.dat" u 1:8 title "7.1", "data_mu_plus_position_8.1.dat" u 1:8 title "8.1", "data_mu_plus_position_8.7.dat" u 1:8 title "8.7", "data_mu_plus_position_9.0.dat" u 1:8 title "9.0", "data_mu_plus_position_9.4.dat" u 1:8 title "9.4", "data_mu_plus_position_9.7.dat" u 1:8 title "9.7", "data_mu_plus_position_9.9.dat" u 1:8 title "9.9"


unset multiplot
