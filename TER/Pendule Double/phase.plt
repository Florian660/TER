set title "Portrait de phase du pendule double"
set autoscale
set xlabel "Position angulaire"
set ylabel "Vitesse angulaire"

u="./euler_explicite.dat"
v="./euler_explicite.dat"

set grid
plot u using 2:3 w l, v using 4:5 w l 
