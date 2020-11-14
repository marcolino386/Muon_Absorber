#plot data about momentum in 4 windows

#set terminal wxt size 1000,900 enhanced font 'Verdana,10' persist

set terminal postscript enhanced color solid  "Times" 16
set output 'figamr0.ps'

set multiplot layout 2,2 rowsfirst


#gráfico a
set label 1
set key top left
set xlabel 'Energy(GeV)'
set ylabel 'P_t {/Symbol D}{P_t}  (GeV/c)'
set xrange[-1:11]
set yrange[-.5:2]
plot "data_mu_plus1_momentum.dat" u 1:4 title "Pt_{/Symbol m+}(0)", "data_mu_minus1_momentum.dat" u 1:4  title "Pt_{/Symbol m-} (0)","data_mu_plus1_momentum.dat" u 1:8 title "Pt_{/Symbol m+}", "data_mu_minus1_momentum.dat" u 1:8 title "Pt_{/Symbol m-}", "data_mu_minus1_momentum.dat" u 1:12 title "{/Symbol D}_{Pt}-", "data_mu_plus1_momentum.dat" u 1:12 title "{/Symbol D}_{Pt}+"

#grafico 2

set label 1 
set xlabel 'Energy(GeV)'
set ylabel 'pt med' 
plot "data_mu_plus1_momentum.dat" u 1:12 title "mu+", "data_mu_minus1_momentum.dat" u 1:12  title "mu-"

#grafico 3

set label 1 
set xlabel 'Energy(GeV)'
set ylabel '{/Symbol D}P_{T}' 
plot "data_mu_plus1_momentum.dat" u 1:12 title "mu+", "data_mu_minus1_momentum.dat" u 1:12  title "mu-"

#gráfico 4
set yrange [0:35000]
set label 1 
set xlabel 'Energy(GeV)'
set ylabel 'N'
plot "data_mu_plus1_momentum.dat" u 1:18 title "mu+", "data_mu_minus1_momentum.dat" u 1:18  title "mu-"

unset multiplot
