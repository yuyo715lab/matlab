set terminal postscript eps enhanced color
set output "pvcurve-N8-big0001.eps"
set size 0.7,0.7
set xrange [4.65:4.675]
set yrange [0.6:0.66]
set grid
set xlabel "Ps[p.u.]"
set ylabel "Vr[p.u.]"
plot "../csv/Node8-01.csv" using 1:2 t 'Node8'
  