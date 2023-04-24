set title "Courbe de convergence du schema Euler explicite"
set xlabel "Pas de temps"
set ylabel "Erreur globale"



mc(x)=a*x+b 
fit mc(x) 'erreur.dat' using (log($1)):(log($2)) via a,b 
plot "erreur.dat" using 1:2, exp(mc(log(x))) title sprintf("droite pente ordre %.3f", a)




