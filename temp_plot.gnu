set autoscale # scale axes automatically
 unset log # remove any log-scaling
 unset label # remove any previous labels
 set xtic auto # set xtics automatically
 set ytic auto # set ytics automatically
 set title "Temperature of the ICM"
 set xlabel "r(kpc)"
 set ylabel "T(K)"
 set logscale x
 set logscale y
 set key at 0.01,100
 #set label "Yield Point" at 0.003,260
 # set arrow from 0.0028,250 to 0.003,280
 #set xr [-2.0:2.0]
 #set yr [-2.0:2.0]
 set key outside

 plot "temperature.dat" using 1:2 title 'ICM Temperature' with lines lt rgb "black"

