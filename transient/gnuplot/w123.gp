set terminal postscript eps enhanced color
set output "../../../../Paper/meeting/transient-analysis/figure/w123.eps"
set grid
set xlabel "time [sec]"
set ylabel "frequency [Hz]"
set size 0.8,0.8

plot "../csv/w123.csv" u 1:2 w l lw 6 t 'G1',\
     "../csv/w123.csv" u 1:3 w l lw 6 t 'G2',\
     "../csv/w123.csv" u 1:4 w l lw 6 t 'G2'