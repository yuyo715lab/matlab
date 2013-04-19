set terminal postscript eps enhanced color
set output "capa-big.eps"
set xrange [4:5]
set yrange [0.35:1]
set key outside
set grid
set xlabel "Ps[p.u.]"
set ylabel "Vr[p.u.]"
plot "../csv/capa-0.1.csv" u 1:2 t '0.1[p.u.]',\
     "../csv/capa-0.3.csv" u 1:2 t '0.3[p.u.]',\
     "../csv/capa-0.5.csv" u 1:2 t '0.5[p.u.]',\
     "../csv/capa-0.8.csv" u 1:2 t '0.8[p.u.]',\
     "../csv/capa-1.1.csv" u 1:2 t '1.1[p.u.]'