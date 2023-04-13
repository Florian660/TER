program resolution
    use constantes_mod
    use edo_mod
    use schema_mod
    implicit none

    type(edo_type)                      :: edo
    type(rk_type)                       :: rk
    real(PR), dimension(:), allocatable :: u,fe
    integer                             :: i,n
    real(PR)                            :: th, vit, delt,xr,yr,x,y,s,E

    ! call init_edo(edo)
    !
    ! t_0 = 0._PR
    !
    ! allocate(u(edo%dim))
    ! allocate(fe(edo%dim))
    ! !test ca marche
    ! !u = (/2._PR,3._PR/)
    ! u = edo%u0
    !
    ! fe = f(u)
    ! !print*, fe
    ! call init_rk(rk,edo) !ca marchepour tous sch
    ! !call printMat(rk%A)
    ! !print*, rk%b
    !
    ! call un_pas_temps(u,edo,rk) !fonctionne correctement
    ! !print*, u


    ! RK 2 ou 4 à choisir dans parametresrk 4e lignes

    !calcul par RK (2ou4) de theta, la vitesse et les coordo cartésiennes
    open(unit=10,file='RK.dat')
    call init_edo(edo)
    allocate(u(2))
    u = edo%u0
    call init_rk(rk,edo)
    delt = edo%tfin/edo%n
        do i = 0,edo%n
            x = l*sin(u(2))
            y = -l*cos(u(2))
            write(10,*) i*delt, u(2), u(1),x,y !temps, th, vitesse,x,y
            call un_pas_temps(u,edo,rk,delt)
        end do
    close(10)

    !cas réel
    ! open(unit = 11, file ='cas_reel.dat')
    ! delt = edo%tfin/edo%n
    ! do i = 0,edo%n
    !     th = threel(edo%u0(2),delt*i)
    !     vit = vitreel(edo%u0(2), delt*i)
    !     xr=l*sin(th)
    !     yr = -l*cos(th)
    !     write(11,*) i*delt, th, vit,xr,yr
    ! end do
    ! close(11)
    !
    !
    ! open(unit = 12, file = 'erreur.dat')
    ! n = 800
    ! do while (n<=10000)
    !     delt = edo%tfin/n
    !     u = edo%u0
    !     th = 0._PR
    !     E= 0._PR
    !     s = 0._PR
    !
    !     do i = 0,n
    !         th = threel(edo%u0(2),i*delt)
    !         s = s + (th-u(2))**2
    !         call un_pas_temps(u,edo,rk,delt) ! pb le deltat est direct dedans faut voir pour le changer
    !
    !
    !     end do
    !
    !     E = sqrt(s*delt)
    !     write(12,*) delt,E
    !     print*, E
    !
    !     n = n*2
    ! end do
    ! close(12)




    ! deallocate(fe)
    deallocate(u)
    call free_edo(edo)
    call free_rk(rk)


end program resolution
