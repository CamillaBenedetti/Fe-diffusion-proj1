set autoscale # scale axes automatically
 unset log # remove any log-scaling
 unset label # remove any previous labels
 set xtic auto # set xtics automatically
 set ytic auto # set ytics automatically
 set title "Fe abundance of the gas cluster at 5 Gyr"
 set xlabel "r (kpc)"
 set ylabel "Z_Fe/Z_S_u_n)"
 set logscale x
 #set logscale y
 set key at 0.01,100
 set key outside

 plot "zfe_initial.dat" using 1:2 title 'Fe abundance at time0' with lines lt rgb "black"
 replot "zfe_initial.dat" using 1:3 title 'Fe abundance obserrved -zout' with lines lt rgb "blue"
 replot "Fe_5Gyr.dat" using 1:2 title 'Fe abundance at 5 Gyr' with lines lt rgb "purple"
 replot "zfe_initial.dat" using 1:7 title 'Fe abundance observed at 5 Gyr' with lines lt rgb "red"
 

