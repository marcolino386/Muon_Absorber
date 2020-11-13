#set terminal wxt size 1000,900 enhanced font 'Verdana,10' persist

set terminal postscript enhanced color solid  "Times" 14
set output 'xyE20kpB00.eps'

set multiplot layout 2,2 rowsfirst


#gráfico a
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

#grafico 3
set label 1
set key top left
set xlabel 'Energy(GeV)'
set ylabel 'E (GeV)'
set xrange[-1:11]
set yrange[-.1:10]
plot  "data_mu_plus2_position.dat" u 1:1 title "E","data_mu_plus2_position.dat" u 1:6 title "<E>", "data_mu_plus2_position.dat" u 1:9  title "{/Symbol D}E"


#gráfico 4
set label 1
set key top left
set xlabel 'Energy(GeV)'
set ylabel '{/Symbol s}'
set xrange[-1:11]
set yrange[-.1:.6]
plot "data_mu_plus2_position.dat" u 1:10 title "{/Symbol s}_x", "data_mu_plus2_position.dat" u 1:11  title "{/Symbol s}_y", "data_mu_plus2_position.dat" u 1:12 title "{/Symbol s}_E"

unset multiplot
