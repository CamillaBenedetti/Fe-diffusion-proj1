reset
set autoscale # scale axes automatically
 unset log # remove any log-scaling
 unset label # remove any previous labels
 set xtic auto # set xtics automatically
 set ytic auto # set ytics automatically
 set title "Fe density in the gas cluster during the evolution"
 set xlabel "r (kpc)"
 set ylabel "rhoFe (g/cm^3)"
 set logscale x
 set logscale y
 set key at 0.01,100
 set key outside

 plot "Fe_1Gyr.dat" using 1:3 title 'rhoFe at 1 Gyr' with lines lt rgb "red"
 replot "Fe_2Gyr.dat" using 1:3 title 'rhoFe mass at 2 Gyr' with lines lt rgb "purple"
 replot "Fe_5Gyr.dat" using 1:3 title 'rhoFe mass at 5 Gyr' with lines lt rgb "blue"
 

