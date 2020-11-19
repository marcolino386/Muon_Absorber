reset
set terminal postscript enhanced color solid  "Times" 10
set output 'ptthN20kpB00.eps'

set multiplot layout 2,2 rowsfirst

#gráfico a
set label 1
set key top left
set xlabel 'Energy(GeV)'
set ylabel 'P_t {/Symbol D}{P_t}  (GeV/c)'
set xrange[-1:11]
set yrange[-.5:1.5]
plot "data_mu_plus2_momentum.dat" u 1:5 title "Pt_{/Symbol m+}(0)",  "data_mu_minus2_momentum.dat" u 1:5  title "Pt_{/Symbol m-} (0)","data_mu_plus2_momentum.dat" u 1:8 title "Pt_{/Symbol m+}", "data_mu_minus2_momentum.dat" u 1:8 title "Pt_{/Symbol m-}", "data_mu_plus2_momentum.dat" u 1:12 title "{/Symbol D}_{Pt}-", "data_mu_minus2_momentum.dat" u 1:12 title "{/Symbol D}_{Pt}+"

#grafico 2
set key top right
set label 1 
set xlabel 'Energy(GeV)'
set ylabel '{/Symbol q} {/Symbol D}{{/Symbol q}}  (grad)'
set xrange[-1:11]
set yrange[-1:15]
plot "data_mu_plus2_momentum.dat" u 1:4 title "{/Symbol q}_{/Symbol m+}(0)", "data_mu_minus2_momentum.dat" u 1:4  title "{/Symbol q}_{/Symbol m-} (0)", "data_mu_plus2_momentum.dat" u 1:9 title "{/Symbol q}_{/Symbol m+}", "data_mu_minus2_momentum.dat" u 1:9 title "{/Symbol q}_{/Symbol m-}", "data_mu_minus2_momentum.dat" u 1:13 title "{/Symbol D}_{{/Symbol q}}-", "data_mu_plus2_momentum.dat" u 1:13 title "{/Symbol D}_{{/Symbol q}}+"

#grafico 3

set label 1 
set key top left
set yrange[-1:5]
set xrange[-2:11]
set xlabel 'Energy(GeV)'
set ylabel '{/Symbol s}_{P_{t}} , {/Symbol s}_{/Symbol q}'
plot "data_mu_minus2_momentum.dat"  u ($1):(20*($23))  title "20*{/Symbol s}_{Pt}-","data_mu_plus2_momentum.dat" u ($1):(20*($23))   title "20*{/Symbol s}_{Pt}+", "data_mu_minus2_momentum.dat" u 1:22   title "{/Symbol s}_{{/Symbol q}}-","data_mu_plus2_momentum.dat" u 1:22   title "{/Symbol s}_{{/Symbol q}}+"


#gráfico 4
set key bottom right
set xrange[-2:11]
set yrange [-500:22000]
set label 1 
set xlabel 'Energy(GeV)'
set ylabel 'N'
plot "data_mu_plus2_momentum.dat" u 1:24 title "N_{/Symbol m+}", "data_mu_minus2_momentum.dat" u 1:24  title "N_{/Symbol m-}"

unset multiplot
