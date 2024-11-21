set autoscale # scale axes automatically
 unset log # remove any log-scaling
 unset label # remove any previous labels
 set xtic auto # set xtics automatically
 set ytic auto # set ytics automatically
 set title "Fe abundance of the gas cluster during the evolution"
 set xlabel "r (kpc)"
 set ylabel "Z_Fe/Z_S_u_n)"
 set logscale x
 #set logscale y
 set key at 0.01,100
 set key outside

 plot "zfe_initial.dat" using 1:2 title 'Fe abundance at time0' with lines lt rgb "black"
 replot "Fe_1Gyr.dat" using 1:2 title 'Fe abundance at 1 Gyr' with lines lt rgb "red"
 replot "Fe_2Gyr.dat" using 1:2 title 'Fe abundance at 2 Gyr' with lines lt rgb "purple"
 replot "Fe_5Gyr.dat" using 1:2 title 'Fe abundance at 5 Gyr' with lines lt rgb "blue"
 replot "initial.dat" using 1:3 title 'Observed Fe abundance' with lines lt rgb "gray"
 

 
 

