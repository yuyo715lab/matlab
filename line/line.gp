set terminal postscript eps enhanced color
set output "line.eps"
set xrange [4.4:4.85]
set yrange [0.4:0.85]
#set key outside
set grid
set xlabel "Ps[p.u.]"
set ylabel "Vr[p.u.]"
plot "./line-0.002.csv" u 1:2 t '-0.002[p.u.]',\
     "./line-0.005.csv" u 1:2 t '-0.005[p.u.]',\
     "./line-0.008.csv" u 1:2 t '-0.008[p.u.]'
