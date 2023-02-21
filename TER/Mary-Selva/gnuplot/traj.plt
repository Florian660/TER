set nokey 
set title "Pendule simple - Euler Explicite - y=f(x)"
unset autoscale x
unset autoscale y
set xrange [-0.07:0.07]
set yrange [-0.13:-0.05]
u = "./euler.dat"
set grid
plot u using 4:5 with linespoints
