set autoscale # scale axes automatically
 unset log # remove any log-scaling
 unset label # remove any previous labels
 set xtic auto # set xtics automatically
 set ytic auto # set ytics automatically
 set title "Density distribution of the galaxy cluster"
 set xlabel "r (kpc)"
 set ylabel "rho (g/cm^3)"
 set logscale x
 set logscale y
 set key at 0.01,100
 #set label "Yield Point" at 0.003,260
 # set arrow from 0.0028,250 to 0.003,280
 #set xr [-2.0:2.0]
 #set yr [-2.0:2.0]
 set key outside

 plot "density.dat" using 1:2 title 'Navarro, Frenk, White profile Density with BCG and Temperature gradient' with lines lt rgb "purple"
 replot "density.dat" using 1:3 title 'Dark Matter Density, analytical formula' with lines lt rgb "red"
 replot "density.dat" using 1:4 title 'Navarro, Frenk, White profile Density without BCG, isothermal' with lines lt rgb "blue"
 #replot "mdensity.dat" using 1:4 title 'Analytical Density' with lines lt rgb "black"
 #missing: Hernquist profile, Density with BCG
