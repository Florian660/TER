program pendules

    implicit none

    integer, parameter :: PR=8
    real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
    integer :: N,i
    real(kind=PR) :: to,thetao,T,g,l,vo,R, h
    real(kind=PR), dimension(:), allocatable :: ti,thetai,vi, coord_x, coord_y


    open(unit=1,file="parametres.dat",action="read")
    read(1,*) g,l,to,T,thetao, vo, N

    allocate(ti(N+1),thetai(N+1),vi(N+1),coord_x(N+1),coord_y(N+1))



    call euler(g,l,to,thetao,vo,T,N,ti,thetai,vi)
    call traject_cart1(l,thetai,coord_x,coord_y)


    open(unit=2,file="euler.dat",action="write")

    do i=1,N+1
        write(2,*) ti(i), thetai(i)*180/PI, vi(i), coord_x(i), coord_y(i)
    end do

    open(unit=3,file="erreur.dat",action="write")
    do i=2,30
        h = 1._PR / (10*2**i)
        write(3,*) h, abs(thetao*cos(sqrt(g/l)*ti(i)) - thetai(i))
    end do


    deallocate(ti,thetai,vi)

    close(1)
    close(2)


    contains

    function f (a,b,c) result(y)
        real(kind=PR), intent(in) :: a,b,c

        real(kind=PR) :: y
        y=-(a/b)*c  ! Cas linéaire petites oscillations
        ! y=-(a/b)*sin(c)

    end function f

    subroutine euler (g,l,to,thetao,vo,T,N,ti,thetai,vi)
        integer, intent(in) :: N ! Nombre de sous-intervalles
        real(kind=PR), intent(in) :: g, l, to, thetao, T, vo

        real(kind=PR), dimension(N+1), intent(out) :: ti,thetai,vi

        integer :: i
        real(kind=PR) :: h



        h=(T-to)/N
        thetai(1)=thetao
        vi(1)=vo
        ti(1)=to

        do i=2,N+1
            ti(i)=ti(i-1)+h
            thetai(i)=thetai(i-1)+h*vi(i-1)
            vi(i)=vi(i-1)+h*f(g,l,thetai(i-1))
        end do

    end subroutine euler

    subroutine traject_cart1(R, theta, coord_x1, coord_y1)
        ! Trace la trajectoire de l'objet dans un repère cartésien pendule simple (initialement polaire).
        real(kind=PR), intent(in) :: R
        real(kind=PR), dimension(N+1), intent(in) :: theta
        real(kind=PR), dimension(N+1), intent(out) :: coord_x1, coord_y1

        coord_x1 = R*sin(theta)
        coord_y1 = -R*cos(theta)

    end subroutine traject_cart1

    subroutine traject_cart2(l1, l2, theta, psi, coord_x2, coord_y2)
        ! Trace la trajectoire de la deuxième masse du pendule double dans un repère cartésien (initialement polaire).
        real(kind=PR), intent(in) :: l1, l2
        real(kind=PR), dimension(N+1), intent(in) :: theta, psi
        real(kind=PR), dimension(N+1), intent(out) :: coord_x2, coord_y2

        coord_x2 = l1*sin(theta)+l2*sin(psi)
        coord_y2 = -(l1*cos(theta)+l2*cos(psi))

    end subroutine traject_cart2


end program
