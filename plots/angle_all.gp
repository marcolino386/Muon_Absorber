set terminal postscript enhanced color solid  "Times" 10
set output 'freq_all.eps'

set xlabel 'Energy(GeV)'
set ylabel 'Frequency'
set xrange[0:10]
set yrange[0:80]
set multiplot

plot "histogram_Energy2_Angle_all_plus.dat" u 1:2 title "2 GeV", "histogram_Energy3_Angle_all_plus.dat" u 1:2 title "3 GeV", "histogram_Energy4_Angle_all_plus.dat" u 1:2 title "4 GeV", "histogram_Energy5_Angle_all_plus.dat" u 1:2 title "5 GeV", "histogram_Energy6_Angle_all_plus.dat" u 1:2 title "6 GeV", "histogram_Energy7_Angle_all_plus.dat" u 1:2 title "7 GeV", "histogram_Energy8_Angle_all_plus.dat" u 1:2 title "8 GeV", "histogram_Energy9_Angle_all_plus.dat" u 1:2 title "9 GeV", "histogram_Energy10_Angle_all_plus.dat" u 1:2 title "10 GeV"
