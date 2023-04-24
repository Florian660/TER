set autoscale
set nokey 
set title "Energie mecanique en fonction du temps"

set xlabel "Temps"
set ylabel "Energie mecanique"

set xrange[0:5]

u="./realite.dat"
v="./euler_explicite.dat"

set grid
plot u using 1:4 w l , v using 1:4 w l 