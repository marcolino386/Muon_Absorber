#set terminal wxt size 1000,900 enhanced font 'Verdana,10' persist

reset
set terminal postscript enhanced color solid  "Times" 10
set output 'fE3e10B05Ea5.711p.eps'
set multiplot layout 1,2 rowsfirst
#gráfico a
set key top left
#histograma
n=1000 #number of intervals
max=3 #max value
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
set xlabel "E (GeV)"
set ylabel "Frequency"
#count and plot
set parametric
plot '/home/lgp/absorber/data/5k/5.7/5.7m/Energy3_Angle_5.711.dat' u (hist($4,width)):(1.0) smooth freq w l t '{/Symbol m+}' ,'/home/lgp/absorber/data/5k/5.7/5.7p/Energy3_Angle_5.711.dat' u (hist($4,width)):(1.0) smooth freq w l t '{/Symbol m-}', "/home/lgp/absorber/data/5k/1.9/1.9m/Energy3_Angle_1.909.dat" u (hist($4,width)):(1.0) smooth freq w l t '{/Symbol m+}' ,"/home/lgp/absorber/data/5k/1.9/1.9p/Energy3_Angle_1.909.dat" u (hist($4,width)):(1.0) smooth freq w l t '{/Symbol m-}',"/home/lgp/absorber/data/5k/2.2/2.2m/Energy3_Angle_2.291.dat" u (hist($4,width)):(1.0) smooth freq w l t '{/Symbol m+}' ,"/home/lgp/absorber/data/5k/2.2/2.2p/Energy3_Angle_2.291.dat" u (hist($4,width)):(1.0) smooth freq w l t '{/Symbol m-}'


#grafico 2
n=1000 #number of intervals
max=10 #max value
min=0. #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set xrange [min:max]
set yrange [0:150]
#to put an empty boundary around the
#data inside an autoscaled graph.
set title '10 GeV'
set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/5,max
set boxwidth width*0.9
set style fill solid 0.5 #fillstyle
set tics out nomirror
set xlabel "E (GeV)"
set ylabel "Frequency"
#count and plot
plot "/home/lgp/absorber/data/5k/5.7/5.7m/Energy10_Angle_5.711.dat" u (hist($4,width)):(1.0) smooth freq w l t '{/Symbol m+}' ,"/home/lgp/absorber/data/5k/5.7/5.7p/Energy10_Angle_5.711.dat" u (hist($4,width)):(1.0) smooth freq w l t '{/Symbol m-}', "/home/lgp/absorber/data/5k/1.9/1.9m/Energy10_Angle_1.909.dat" u (hist($4,width)):(1.0) smooth freq w l t '{/Symbol m+}' ,"/home/lgp/absorber/data/5k/1.9/1.9p/Energy10_Angle_1.909.dat" u (hist($4,width)):(1.0) smooth freq w l t '{/Symbol m-}', "/home/lgp/absorber/data/5k/2.2/2.2m/Energy10_Angle_2.291.dat" u (hist($4,width)):(1.0) smooth freq w l t '{/Symbol m+}' ,"/home/lgp/absorber/data/5k/2.2/2.2p/Energy10_Angle_2.291.dat" u (hist($4,width)):(1.0) smooth freq w l t '{/Symbol m-}'


unset multiplot
