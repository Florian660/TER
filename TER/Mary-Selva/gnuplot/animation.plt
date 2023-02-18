reset
set nokey       # Désactive la légende
set term gif animate delay 4 size 854,480   # délai de (4*0.01)s entre chaque image => 1/0.04=25 images/sec, gif de 854*480 pixels

set output 'traj_pend_simple.gif'       # nom du fichier .gif

set xr[-40:40]      # xrange
set yr[-40:40]      # yrange
set xl 'x' font 'Times New Roman:Italic, 20'        # xlabel
set yl 'y' font 'Times New Roman:Italic, 20'        # ylabel
set tics font 'Times New Roman,18'                  # graduations

set size square     # graphe carré
set grid

#-----------------------PARAMETRES-----------------------
m = 0.1             # en kg
l  = 20             # en m
g  = 9.81           # en m.s²
theta = 1.0472      # en rad        (Exemples: 10° = 0.1745 rad  ou  60° = 1.0472 rad)
dtheta = 0          # en s^-1
N = 400000.0        # nb de sous-intervalles
t = 0               # t0
tn = 4000.0         # t final (Exemple 4000.0 => environ 10 sec de simulation)
dt = (tn-t)/N       # pas de temps
R  = 1              # rayon de la masse du pendule
r  = 0.1            # rayon des points de la trajectoire

#-----------------------LABELS-----------------------
label(a, b, c, d, e) = sprintf("\
{/Times:Italic m} = %.2f kg\n\
{/Times:Italic l}   = %3.1f m\n\
{/Times:Italic g}    = %3.2f m/s^2\n\
{/Times:Italic dt}   = %.2f s\n\
{/symbol-oblique q_{/Times:Normal 0}}  = %3.2f °\n", a, b, c, d, e)

# sprintf convertit en chaîne de caractères et format "%5.2f" = champ à 5 espaces (2 décimales, 1 pt, 2 chiffres partie entière)

#-----------------------POSITION INITIALE-----------------------
# Première boucle pour figer l'image à la position initiale
do for [i = 1:20] {

    set title "Mouvement pendule simple - Euler Explicite" font 'Times:Normal, 20'

    # Fil
    set arrow 1 nohead lw 2 from 0, 0 to l*sin(theta), -l*cos(theta) lc -1
    # flèche 1 sans tête, lw = line_width 2 fois plus large, vecteur de 0 à (l*sin(theta), -l*cos(theta)), lc -1 = noire

    # Paramètres
    set label 1 left at 45, 20 label(m, l, g, dt, theta*180/pi) font 'Times:Normal, 18'     # label des paramètres à la position (45,20)

    # Liaison pivot
    set object 1 circle at 0, 0 fc rgb 'black' size R fs solid front

    # Masse
    set object 2 circle at l*sin(theta), -l*cos(theta) fc rgb 'blue' size R fs solid front      # fc=fillcolor, fs=fillstyle solid=remplissage avec couleur unie, front=mettre au premier plan

    # Trace
    plot 1/0   # 1/0 n'est pas défini, permet de tracer seulement les objets utiles

    # Trajectoire
    set object 2 size r     # réduit le cercle de la masse suspendue pour créer une trajectoire
}

# Équation caractéristique pendule simple
f(a,b,c) = -(a/b)*sin(c)

#-----------------------ANIMATION-----------------------
do for [i = 1:tn] {
    t = t + dt

    # Euler Explicite
    theta = theta + dt*dtheta
    dtheta = dtheta + dt*f(g,l,theta)

    # Mouvement du fil
    set arrow 1 nohead lw 2 from 0, 0 to l*sin(theta), -l*cos(theta) lc -1

    # Mouvement de la masse
    set object i+2 circle at l*sin(theta), -l*cos(theta) fc rgb "blue" size R fs solid front

    # Trace et retire les points inutiles
    if (i%20==0){
        plot 1/0
    }

    # Trajectoire
    set object i+2 size r     # réduit le cercle de la masse suspendue pour créer une trajectoire
}

set out
