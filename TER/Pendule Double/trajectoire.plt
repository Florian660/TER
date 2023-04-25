set title "Trajectoire du pendule double"
set autoscale 
set xlabel "x"
set ylabel "y"


u ="./euler_explicite.dat"
v ="./euler_explicite.dat"



set grid 
plot u using 6:7 w l, v using 8:9 w l  