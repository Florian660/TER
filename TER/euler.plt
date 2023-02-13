set title "euler"
set autoscale 
set xlabel "time in sec "
set ylabel "x"
m= "./euler.dat"
set grid 
plot m using 1:2 with linespoints 