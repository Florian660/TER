
set title "Trajectoire du pendule simple"
set autoscale 
set xlabel "x"
set ylabel "y"


u ="./euler_explicite.dat"



set grid 
plot u using 5:6 w l