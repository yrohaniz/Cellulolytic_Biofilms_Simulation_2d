set terminal postscript eps color enhanced "Times-Roman" 24
set size 2.2,2.2
set size square
set key out
set xrange [0:24]
set xlabel "time"
set ylabel "Total Carbon" 
set yrange [0.2:1.2]
set xtics 0,1,24
set tics out nomirror
set title "Total amount of carbon substrate for different values of the coefficient of the Wiener process vs. time"
xi = "0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00"
set output 'carbon.eps'

plot for [i in xi] sprintf('Spatially_integrated_carbconcn_512x512_xi-%s', i) w l lw 3 ti 'xi='.i
