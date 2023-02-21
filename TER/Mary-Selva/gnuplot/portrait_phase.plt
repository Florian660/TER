set title "Pendule simple - Euler - dtheta=f(theta)"
set autoscale
set xlabel "Theta (deg)"
set ylabel "Dtheta (s-1)"
u = "./euler.dat"
set grid
plot u using 2:3 with linespoints
