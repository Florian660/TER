set xlabel "t"
set ylabel "theta1, theta2"
plot "pendule_couple.dat" using 1:2 with lines title "theta1", \
     "pendule_couple.dat" using 1:3 with lines title "theta2"
