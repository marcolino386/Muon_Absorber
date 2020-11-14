#set terminal wxt size 1000,900 enhanced font 'Verdana,10' persist

reset
set terminal postscript enhanced color solid  "Times" 10
set output 'freqxy3e10GeV1kpB05.eps'
set multiplot layout 2,2 rowsfirst

# linha
set parametric
const=0.647957
constxm3=0.531577
constym3=0.681741
set trange [0:100]
#gráfico a
set key top left
#histograma
n=1000 #number of intervals
max=1.3 #max value
min=0. #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set xrange [min:max]
set yrange [0:100]
#to put an empty boundary around the
#data inside an autoscaled graph.
set title '3 GeV'
set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/5,max
set boxwidth width*0.9
set style fill solid 0.5 #fillstyle
set tics out nomirror
set xlabel "x (m)"
set ylabel "Frequency"
#count and plot
set parametric
plot "Energy3_Angle_7.014.dat" u (hist($2,width)):(1.0) smooth freq w boxes lc rgb"green" t '' ,"Energy3_Angle_7.014.dat" u (hist($3,width)):(1.0) smooth freq w l t ''


#grafico 2
#reset
# linha
set parametric
const=0.647957
constxm=0.614384
constym=0.681741
set trange [0:350]
#histograma
n=1000 #number of intervals
max=1.3 #max value
min=0. #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set xrange [min:max]
set yrange [0:400]
#to put an empty boundary around the
#data inside an autoscaled graph.
set title '10 GeV'
set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/5,max
set boxwidth width*0.9
set style fill solid 0.5 #fillstyle
set tics out nomirror
set xlabel "x (m)"
set ylabel "Frequency"
#count and plot
plot "Energy10_Angle_7.014.dat" u (hist($2,width)):(1.0) smooth freq w boxes lc rgb"green" t '' ,"Energy10_Angle_7.014.dat" u (hist($3,width)):(1.0) smooth freq w l t ''



#grafico 3
set label 1
set title '3 GeV'
set key top right
set xlabel 'x (m)'
set ylabel 'y (m)'
set xrange [min:max]
set xtics min,(max-min)/5,max
set yrange[0:1.3]
plot "Energy3_Angle_7.014.dat" u 2:3 title "E = 3 GeV"


#gráfico 4
set label 1
set title '10 GeV'
set key top right
set xlabel 'x (m)'
set ylabel 'y (m)'
set xrange [min:max]
set xtics min,(max-min)/5,max
set yrange[0:1.3]
plot "Energy10_Angle_7.014.dat" u 2:3 title "E = 10 GeV"



unset multiplot
