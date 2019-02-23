set terminal postscript eps color enhanced "Times-Roman" 24
set size 2.2,2.2
set size square
set xrange [0:1]
set xlabel "X"
set ylabel "Y" 
set yrange [0:1]
set cbrange [0:1]
set xtics 0,0.1,1
set ytics 0,0.1,1
set tics out nomirror

t = 60
n = "512x512"

set output sprintf("%s_t".t.".eps", n)
set title 'Biomass load on the substrate captured at t='.t
plot sprintf('Biomass_density_%s_formatted_t'.t, n) matrix nonuniform with image notitle


