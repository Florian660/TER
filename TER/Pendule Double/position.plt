set title "Position angulaire du pendule double en fonction du temps"
set autoscale 
set xlabel "Temps"
set ylabel "Position angulaire"

m="./euler_explicite.dat"
n="./euler_explicite.dat"


set grid 
plot m using 1:2 w l, n using 1:4 w l