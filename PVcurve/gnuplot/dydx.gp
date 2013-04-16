#set terminal postscript eps enhanced color
#set output "dvdp-test.eps"
set key top right
set size 0.7,0.7
set log y
set grid
set yrange [100:0.01]
set xlabel 'Ps [p.u.]'
set ylabel '|dV/dP|'
plot "../csv/dvdp2.csv" u 1:(-1*($2)) t 'dV/dP'