set terminal postscript eps enhanced color
set output "../../../../Paper/meeting/transient-analysis/figure/v123.eps"
set grid
set xlabel "Time [sec]"
set ylabel "Voltage [p.u.]"
set size 0.8,0.8

plot "../csv/v123.csv" u 1:2 w l lw 6 t 'G1',\
     "../csv/v123.csv" u 1:3 w l lw 6 t 'G2',\
     "../csv/v123.csv" u 1:4 w l lw 6 t 'G2'