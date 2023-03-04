set nokey 
set title "Pendule simple - Euler - theta=f(t)"
set autoscale
set xlabel "Temps (s)"
set ylabel "Theta (deg)"
u = "./euler.dat"
set grid
theta(x) = 10*cos(sqrt(9.81/0.1)*x)
plot theta(x), u using 1:2 with linespoints
