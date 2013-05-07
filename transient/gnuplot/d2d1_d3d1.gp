set terminal postscript eps enhanced color
set output "../../../../Paper/meeting/transient-analysis/figure/d2d1_d3d1.eps"
set grid
set xlabel "time [sec]"
set ylabel "phase difference angle [degree]"
set size 0.8,0.8

plot "../csv/d2d1_d3d1.csv" u 1:($2)/pi*180 w l lw 6 t '{/Symbol d}2-{/Symbol d}1',\
     "../csv/d2d1_d3d1.csv" u 1:($3)/pi*180 w l lw 6 t '{/Symbol d}3-{/Symbol d}1'