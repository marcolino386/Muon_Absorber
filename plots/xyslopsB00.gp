#set terminal wxt size 1000,900 enhanced font 'Verdana,10' persist

reset
set terminal postscript enhanced color solid  "Times" 10
set output 'xyslopsB00.eps'
set multiplot layout 2,2 rowsfirst


#gráfico a
set label 1
set key top left
set xlabel 'Energy(GeV)'
set ylabel 'SlopeX'
set xrange[-1:11]
set yrange[-.1:.2]
plot "data_mu_plus2_momentum.dat" u 1:16 title "{/Symbol a_2^+}(0)", "data_mu_plus2_momentum.dat" u 1:18  title "{/Symbol {<a_2^+>}}", "data_mu_plus2_momentum.dat" u 1:20 title "{/Symbol s_{a_2^+}}"

#grafico 2
set key top right
set label 1 
set ylabel 'SlopeX'
set xrange[-1:11]
set yrange[-.1:.2]
plot "data_mu_minus2_momentum.dat" u 1:16 title "{/Symbol a_2^-}(0)", "data_mu_minus2_momentum.dat" u 1:18  title "{/Symbol {<a_2^->}}","data_mu_minus2_momentum.dat" u 1:20 title "{/Symbol s_{a_2^-}}"

#grafico 3
set label 1
set key top right
set xlabel 'Energy(GeV)'
set ylabel 'SlopeY'
set xrange[-1:11]
set yrange[-.1:.15]
plot "data_mu_plus2_momentum.dat" u 1:17 title "{/Symbol a_4^+}(0)", "data_mu_plus2_momentum.dat" u 1:19  title "{/Symbol <{a_4^+}>}", "data_mu_plus2_momentum.dat" u 1:21 title "{/Symbol s_{a_4^+}}"


#gráfico 4
set label 1
set key top left
set xlabel 'Energy(GeV)'
set ylabel 'SlopeY'
set xrange[-1:11]
set yrange[-.1:.15]
plot "data_mu_minus2_momentum.dat" u 1:17 title "{/Symbol a_4^-}(0)", "data_mu_minus2_momentum.dat" u 1:19  title "{/Symbol {<a_4^->}}",  "data_mu_minus2_momentum.dat" u 1:21 title "{/Symbol s_{a_4^-}}"

unset multiplot
