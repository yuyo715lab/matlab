set terminal postscript eps enhanced color solid
set output "../../../../Paper/meeting/transient-analysis/figure/delta_4145_00001.eps"
set grid
set xlabel "time [sec]"
set ylabel "phase difference angle [degree]"
set size 0.8,0.8
set key below
plot "../csv/delta_44_00001.csv" u 1:2 w l lw 6 t '0.44[sec]',\
     "../csv/delta_43_00001.csv" u 1:2 w l lw 6 t '0.43[sec]',\
     "../csv/delta_45_00001.csv" u 1:2 w l lw 6 t '0.45[sec]',\
     "../csv/delta_4142_00001.csv" u 1:2 w l lw 6 t '0.41[sec]',\
     "../csv/delta_4142_00001.csv" u 1:3 w l lw 6 t '0.42[sec]'
