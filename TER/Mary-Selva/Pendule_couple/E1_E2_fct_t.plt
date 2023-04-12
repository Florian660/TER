set nokey
set xlabel "Temps (s)"
set ylabel "Energie mecanique (J)"
plot "pendule_couple.dat" using 1:8 with lines title "Energie mecanique 1"
replot "pendule_couple.dat" using 1:9 with lines title "Energie mecanique 1"
