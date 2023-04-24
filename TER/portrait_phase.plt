
set title "Portrait de phase du pendule simple"
set autoscale
set xlabel "Position angulaire"
set ylabel "Vitesse angulaire"

v="./euler_explicite.dat"

set grid
plot v using 2:3 w l 

