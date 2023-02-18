J'ai repris le code de Florian et j'ai ajouté une subroutine pour tracer la trajectoire en cartésien. Il y a donc deux nouvelles 
colonnes au fichier euler.dat avec les coordonnées en x et en y.

Aussi, j'ai traité le cas général (sans l'approximation des petits angles donc avec le sinus).
Et là il faut faire attention que sin prend en argument des radians, c'est pour cela que dans le fichier parametres.dat,
j'ai changé la condition initiale en radian (10°=0,175rad) et qu'à des endroits dans le code il y a des conversions.

Ce qui est contenu dans le fichier gnuplot doit être lancé dans gnuplot avec la ligne de commande: load "nom_fichier.plt". Le fichier .gif de ce dossier est ce que l'on obtient en exécutant le fichier "animation.plt".
