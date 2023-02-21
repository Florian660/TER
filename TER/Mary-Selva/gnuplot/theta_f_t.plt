set title "Pendule simple - Euler - theta=f(t)"
set autoscale
set xlabel "Temps (s)"
set ylabel "Theta (deg)"
u = "./euler.dat"
set grid
plot u using 1:2 with linespoints