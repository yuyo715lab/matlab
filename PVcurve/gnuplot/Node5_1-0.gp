set terminal postscript eps enhanced color
set output "pvcurve-N8-big01.eps"
set xrange [4.65:4.675]
set yrange [0.6:0.66]
set grid
set xlabel "Ps[p.u.]"
set ylabel "Vr[p.u.]"
plot "../csv/Node8.csv" u 1:2 t 'Node8'