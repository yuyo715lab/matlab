set terminal postscript eps enhanced color
set output "pvcurve00.eps"

set grid
set xlabel "Ps"
set ylabel "Vr"
plot "../csv/Node5_1-0.csv" using 1:2 t 'Node5',\
     "../csv/Node6.csv" u 1:2 t 'Node6',\
     "../csv/Node8.csv" u 1:2 t 'Node8'