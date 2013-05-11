#set terminal postscript eps enhanced color solid
#set output "../../../../Paper/meeting/transient-analysis/figure/delta_eftime_detail.eps"
set grid
set xlabel "time [sec]"
set ylabel "phase difference angle [degree]"
set size 0.8,0.8
set key below
plot "../csv/test.csv" u 1:2 w l lw 6 t '0.2[sec]',\
     "../csv/test.csv" u 1:4 w l lw 6 t '0.4[sec]',\
     "../csv/test.csv" u 1:5 w l lw 6 t '0.5[sec]'
#          "../test.csv" u 1:3 w l lw 6 t '0.3[sec]',\
#     "../test.csv" u 1:6 w l lw 6 t '0.6[sec]'
#     "../csv/test.csv" u 1:7 w l lw 6 t '0.45[sec]',\
#     "../csv/test.csv" u 1:8 w l lw 6 t '0.46[sec]',\
#     "../csv/test.csv" u 1:10 w l lw 6 t '0.48[sec]',\
#     "../csv/test.csv" u 1:10 w l lw 6 t '0.50[sec]'