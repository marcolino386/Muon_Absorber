#set terminal wxt size 1000,900 enhanced font 'Verdana,10' persist

set terminal postscript enhanced color solid  "Times" 14
set output 'xyE1kpB05.eps'

set multiplot layout 2, rowsfirst


#gráfico 1
set label 1
set key top left
set xlabel 'Energy(GeV)'
set ylabel 'x(m)'
set xrange[-1:11]
set yrange[-.1:1]
plot "data_mu_plus2_position.dat" u 1:2 title "x(0)", "data_mu_plus2_position.dat" u 1:4  title "<x>", "data_mu_plus2_position.dat" u 1:7 title "{/Symbol D}x"

#grafico 2
set key top right
set label 1 
set ylabel 'y(m)'
set xrange[-1:11]
set yrange[-.1:1]
plot "data_mu_plus2_position.dat" u 1:3 title "y(0)", "data_mu_plus2_position.dat" u 1:5  title "<y>", "data_mu_plus2_position.dat" u 1:8 title "{/Symbol D}y"


unset multiplot
