program pendule_couple

    implicit none

    integer, parameter :: n = 50000 ! nombre d'iterations
    integer, parameter :: PR = 8
    real(PR), parameter :: pi = 4._PR*atan(1._PR)
    integer :: i
    real(PR) :: t, dt
    real(PR) :: theta1, vit1, acc1
    real(PR) :: theta2, vit2, acc2
    real(PR) :: L1, L2, m1, m2, C, g
    real(PR) :: E, E1, E2, Ec1, Ec2, Ep1, Ep2
    real(PR) :: x1, y1, x2, y2

    ! Parametres
    L1 = 1._PR  ! longueur du premier fil
    L2 = 1._PR  ! longueur du deuxieme fil
    m1 = 1._PR  ! masse de la premiere masse
    m2 = 1._PR  ! masse de la deuxieme masse
    C = 1._PR   ! constante de couplage
    g = 9.81    ! acceleration de la gravite

    ! Conditions initiales / theta en rad
    theta1 = 10._PR * pi/180    ! Angle 1
    vit1 = 0._PR                ! Vitesse 1
    acc1 = 0._PR                ! Accélération 1
    theta2 = 0._PR * pi/180     ! Angle 2
    vit2 = 0._PR                ! Vitesse 2
    acc2 = 0._PR                ! Accélération 2

    ! Pas de temps
    dt = 0.01

    ! Ouverture du fichier de sortie
    open(1, file= 'pendule_couple.dat', status="replace")

    ! Euler explicite
    do i = 1, n

        t = i * dt

        ! Calcul de alpha1 et alpha2
        acc1 = (-g/L1) * sin(theta1) + (C/(m1*L1**2)) * (theta2 - theta1)
        acc2 = (-g/L2) * sin(theta2) + (C/(m2*L2**2)) * (theta1 - theta2)

        ! Mise a jour de omega1 et omega2
        vit1 = vit1 + acc1 * dt
        vit2 = vit2 + acc2 * dt

        ! Mise a jour de theta1 et theta2
        theta1 = theta1 + vit1 * dt
        theta2 = theta2 + vit2 * dt

        ! Conversion des coordonnees polaires en coordonnees cartesiennes
        x1 = L1 * sin(theta1)
        y1 = -L1 * cos(theta1)
        x2 = L2 * sin(theta2)
        y2 = -L2 * cos(theta2)

        ! Calcul de l'énergie cinétique et de l'énergie potentielle
        Ec1 = 0.5 * m1 * vit1**2
        Ec2 = 0.5 * m2 * vit2**2
        Ep1 = - m1 * g * L1 * cos(theta1)
        Ep2 = - m2 * g * L2 * cos(theta2)

        ! Calcul des énergies
        E1 = Ec1 + Ep1
        E2 = Ec2 + Ep2

        ! Ecriture des donnees dans le fichier
        write(1, *) t, theta1, theta2, x1, y1, x2, y2, E1, E2, E

    end do

    ! Fermeture du fichier de sortie
    close(1)

end program pendule_couple
