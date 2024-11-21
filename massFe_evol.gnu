set autoscale # scale axes automatically
 unset log # remove any log-scaling
 unset label # remove any previous labels
 set xtic auto # set xtics automatically
 set ytic auto # set ytics automatically
 set title "Fe mass in the gas cluster during the diffusion evolution"
 set xlabel "r(kpc)"
 set ylabel "M_F_e/M_S_u_n"
 set logscale x
 #set logscale y #occhio perché per z=0 non riesce a plottare log(0)
 set key at 0.01,100
 set key outside

 plot "zfe_initial.dat" using 1:5 title 'initial Fe mass' with lines lt rgb "black"
 replot "Fe_1Gyr.dat" using 1:4 title 'Fe mass at 1 Gyr' with lines lt rgb "red"
 replot "Fe_2Gyr.dat" using 1:4 title 'Fe mass at 2 Gyr' with lines lt rgb "purple"
 replot "Fe_5Gyr.dat" using 1:4 title 'Fe mass at 5 Gyr' with lines lt rgb "blue"
 