set autoscale # scale axes automatically
 unset log # remove any log-scaling
 unset label # remove any previous labels
 set xtic auto # set xtics automatically
 set ytic auto # set ytics automatically
 set title "Mass distribution of the galaxy cluster"
 set xlabel "r(kpc)"
 set ylabel "M/M_S_u_n"
 set logscale x
 set logscale y
 set key at 0.01,100
 #set label "Yield Point" at 0.003,260
 # set arrow from 0.0028,250 to 0.003,280
 #set xr [-2.0:2.0]
 #set yr [-2.0:2.0]
 set key outside

 plot "mass.dat" using 1:2 title 'Navarro, Frenk, White profile Mass' with lines lt rgb "purple"
 replot "mass.dat" using 1:3 title 'Dark Matter Mass, analytical formula' with lines lt rgb "red" 
 replot "mass.dat" using 1:4 title 'Hernquist profile, Luminous Mass with BCG' with lines lt rgb "black"

