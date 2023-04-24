set title "Position angulaire du pendule simple en fonction du temps"
set autoscale 
set xlabel "Temps"
set ylabel "Position angulaire"


n="./euler_explicite.dat"


set grid 
plot n using 1:2 w l
