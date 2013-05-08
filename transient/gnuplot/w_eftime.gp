set terminal postscript eps enhanced color
set output "../../../../Paper/meeting/transient-analysis/figure/w_eftime.eps"
set grid
set xlabel "time [sec]"
set ylabel "frequency [Hz]"
set size 0.8,0.8

plot "../csv/w_eftime.csv" u 1:2 w l lw 6 t '0.2[sec]',\
     "../csv/w_eftime.csv" u 1:3 w l lw 6 t '0.3[sec]',\
     "../csv/w_eftime.csv" u 1:4 w l lw 6 t '0.4[sec]',\
     "../csv/w_eftime.csv" u 1:5 w l lw 6 t '0.5[sec]',\
     "../csv/w_eftime.csv" u 1:6 w l lw 6 t '0.6[sec]'