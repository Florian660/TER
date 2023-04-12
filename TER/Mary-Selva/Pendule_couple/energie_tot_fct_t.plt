set nokey
set xlabel "Temps (s)"
set ylabel "Energie totale (J)"
plot "pendule_couple.dat" using 1:10 with lines title "Energie totale"
