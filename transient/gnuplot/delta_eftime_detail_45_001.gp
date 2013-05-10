set terminal postscript eps enhanced color solid
set output "../../../../Paper/meeting/transient-analysis/figure/delta_eftime_detail_45_0001.eps"
set grid
set xlabel "time [sec]"
set ylabel "phase difference angle [degree]"
set size 0.8,0.8
set key below
plot "../csv/delta_eftime_detail_45_0001.csv" u 1:2 w l lw 6 t '0.40[sec]',\
     "../csv/delta_eftime_detail_45_0001.csv" u 1:3 w l lw 6 t '0.41[sec]',\
     "../csv/delta_eftime_detail_45_0001.csv" u 1:4 w l lw 6 t '0.42[sec]',\
     "../csv/delta_eftime_detail_45_0001.csv" u 1:5 w l lw 6 t '0.43[sec]',\
     "../csv/delta_eftime_detail_45_0001.csv" u 1:6 w l lw 6 t '0.44[sec]',\
     "../csv/delta_eftime_detail_45_0001.csv" u 1:7 w l lw 6 t '0.45[sec]',\
     "../csv/delta_eftime_detail_45_0001.csv" u 1:8 w l lw 6 t '0.46[sec]',\
     "../csv/delta_eftime_detail_45_0001.csv" u 1:10 w l lw 6 t '0.48[sec]',\
     "../csv/delta_eftime_detail_45_0001.csv" u 1:10 w l lw 6 t '0.50[sec]'