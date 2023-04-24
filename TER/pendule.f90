program pendules

    implicit none 

    integer, parameter :: PR=8
    real(PR), parameter :: PI=4*ATAN(1._PR)
    integer :: N,i,j,meth,systeme 
    real(kind=PR) :: to,thetao,T,g,l,m,vo,h,S
    real(kind=PR), dimension(:), allocatable :: ti,thetai,vi,coord_x, coord_y, thetai_c, vi_c,coord_x_c,coord_y_c, err
    real(kind=PR), dimension(:), allocatable :: emi,emic

    !l'argument dans un cosinus ou sinus est de base en radian

    open(unit=1,file="parametres.dat",action="read")
    read(1,*) m,g,l,to,T,thetao,vo,N,meth,systeme 
    !systeme : 1=pendule simple linéarisé ; 2=pendule simple 

    select case(meth) !meth : sol exacte ou schémas numériques 
    case(1) !Solution exacte 
        if (systeme==1) then 
            allocate(ti(N+1),thetai(N+1),vi(N+1),emi(N+1),coord_x(N+1),coord_y(N+1))
            call realite(m,g,l,to,T,thetao,vo,N,ti,thetai,vi,emi,systeme)
            call traject_cart(l,thetai,coord_x,coord_y)
            open(unit=2,file='realite.dat',action='write')
            do i=1,N+1 
                write(2,*) ti(i),thetai(i),vi(i),emi(i),coord_x(i),coord_y(i)
            end do 
        else 
            print*, "on connait que la solution exacte du pendule simple linearise"
        end if 
    case(2) !euler explicite 
        allocate(ti(N+1),thetai(N+1),vi(N+1),emi(N+1),coord_x(N+1),coord_y(N+1))
        call euler_explicite(m,g,l,to,T,thetao,vo,N,ti,thetai,vi,emi,systeme)
        call traject_cart(l,thetai,coord_x,coord_y)
        open(unit=3,file='euler_explicite.dat',action='write')
        do i=1,N+1 
            write(3,*) ti(i),thetai(i),vi(i),emi(i),coord_x(i),coord_y(i)
        end do
        deallocate(ti,thetai,vi,emi,coord_x,coord_y)
    case(3) !euler implicite 
        allocate(ti(N+1),thetai(N+1),vi(N+1),emi(N+1),coord_x(N+1),coord_y(N+1))
        call euler_implicite(m,g,l,to,T,thetao,vo,N,ti,thetai,vi,emi,systeme)
        call traject_cart(l,thetai,coord_x,coord_y)
        open(unit=4,file='euler_implicite.dat',action='write')
        do i=1,N+1 
            write(4,*) ti(i),thetai(i),vi(i),emi(i),coord_x(i),coord_y(i)
        end do
    case default 
        print*, "tapez 1"
    end select 

       
    !CALCUL DE L'ERREUR 
    ! allocate(err(8))
    ! open(unit=5,file='erreur.dat',action="write")
    ! do i=1,8  ! on fait varier le pas de temps de 10^-3 à 10^-5 
    !     S=0
    !     N=50*2**(i+4) !pas de temps de 2.5*10**(-2) à 2.4*10**(-5)
    !     allocate(ti(N+1),thetai(N+1),vi(N+1),coord_x(N+1),coord_y(N+1),thetai_c(N+1),vi_c(N+1),coord_x_c(N+1), &
    !     coord_y_c(N+1),emi(N+1), emic(N+1))  
    !     call euler_implicite(m,g,l,to,T,thetao,vo,N,ti,thetai_c,vi_c,emic,systeme)  !on calcule la "vraie" valeur de theta pour chaque pas de temps
    !     call euler_explicite(m,g,l,to,T,thetao,vo,N,ti,thetai,vi,emi,systeme) ! de même pour le schéma d'euler 
    !     do j=1,N+1
    !         S=S+(thetai_c(j)-thetai(j))**2
    !     end do 
    !     h=(T-to)/N
    !     err(i)=SQRT(S*h)  ! on calcule l'erreur en norme L2
    !     write(5,*) h, err(i)
    !     deallocate(ti,thetai,vi,coord_x,coord_y,thetai_c,vi_c,coord_x_c,coord_y_c,emi,emic)
    ! end do 
    ! deallocate(err)



    contains

    function em(a,b,c,d,e,systeme) result(y)
        integer :: systeme 
        real(kind=PR), intent(in) :: a,b,c,d,e

        real(kind=PR) :: y 
        if (systeme==1) then !linearise 
            y=(1/2._PR)*a*(b**2)*(c**2)+a*b*d*(e**2)/2._PR
        else !non linearise
            y=(1/2._PR)*a*(b**2)*(c**2)+a*b*d*(1-COS(e*PI/180._PR))
        end if 

    end function em

    subroutine euler_explicite (m,g,l,to,T,thetao,vo,N,ti,thetai,vi,emi,systeme)
        integer, intent(in) :: N,systeme !nombre de sous-intervalles 
        real(kind=PR), intent(in) :: m,g, l, to, thetao, T, vo 
        real(kind=PR), dimension(N+1), intent(out) :: ti,thetai,vi,emi

        integer :: i 
        real(kind=PR) :: h

        h=(T-to)/N 
        thetai(1)=thetao
        vi(1)=vo
        ti(1)=to
        emi(1)=em(m,l,vo,g,thetao,systeme)

        if (systeme==1) then !linearise 
            do i=2,N+1
                ti(i)=ti(i-1)+h
                thetai(i)=thetai(i-1)+h*vi(i-1)
                vi(i)=vi(i-1)-h*(g/l)*thetai(i-1)
                emi(i)=em(m,l,vi(i),g,thetai(i),systeme)
            end do 
        else !non linearise
            do i=2,N+1
                ti(i)=ti(i-1)+h
                thetai(i)=thetai(i-1)+h*vi(i-1)
                vi(i)=vi(i-1)-h*(g/l)*SIN(thetai(i-1)*PI/180._PR)
                emi(i)=em(m,l,vi(i),g,thetai(i),systeme)
            end do 
        end if 

    end subroutine euler_explicite 

    subroutine euler_implicite (m,g,l,to,T,thetao,vo,N,ti,thetai,vi,emi,systeme)   !petites oscillations
        integer, intent(in) :: N,systeme !nombre de sous-intervalles 
        real(kind=PR), intent(in) :: m,g,l,to, thetao,T,vo 
        real(kind=PR), dimension(N+1), intent(out) :: ti,thetai,vi,emi

        integer :: i 
        real(kind=PR) :: h,p,r

        p=l+g*(h**2)
        r=g*(h**2)

        h=(T-to)/N 
        thetai(1)=thetao
        vi(1)=vo
        ti(1)=to
        emi(1)=em(m,l,vo,g,thetao,systeme)


        if (systeme==1) then !linearise 
            do i=2,N+1
                ti(i)=ti(i-1)+h
                thetai(i)=(1._PR-(r/p))*thetai(i-1)+(h*l/p)*vi(i-1)
                vi(i)=(l/p)*(vi(i-1)-h*(g/l)*thetai(i-1))
                emi(i)=em(m,l,vi(i),g,thetai(i),systeme)
            end do 
        else !non linearise
            ! do i=2,N+1
            !     ti(i)=ti(i-1)+h
            !     thetai(i)=thetai(i-1)+h*vi(i-1)
            !     vi(i)=vi(i-1)+h*f(g,l,thetai(i-1),systeme)
            !     emi(i)=em(m,l,vi(i),g,thetai(i),systeme)
            ! end do 
        end if 
    end subroutine euler_implicite

    subroutine realite(m,g,l,to,T,thetao,vo,N,ti,thetai,vi,emi,systeme)
        integer, intent(in) :: N,systeme !nombre de sous-intervalles 
        real(kind=PR), intent(in) :: m,g,l,to,T,thetao,vo
        real(kind=PR), dimension(N+1), intent(out) :: ti,thetai,vi,emi

        integer :: i 
        real(kind=PR) :: h,omega

        h=(T-to)/N 
        thetai(1)=thetao
        ti(1)=to
        vi(1)=vo
        emi(1)=em(m,l,vo,g,thetao,systeme)
       
        omega=SQRT(g/l)

        do i=2,N+1
            ti(i)=ti(i-1)+h
            thetai(i)=thetao*COS(omega*ti(i))+(vo/omega)*SIN(omega*ti(i))
            vi(i)=-thetao*omega*SIN(omega*ti(i))+vo*COS(omega*ti(i))
            emi(i)=em(m,l,vi(i),g,thetai(i),systeme)
        end do

    end subroutine realite


    subroutine traject_cart(R,theta,coord_x,coord_y)
        ! Trace la trajectoire de l'objet dans un repère cartésien (initialement polaire).
        real(kind=PR), intent(in) :: R
        real(kind=PR), dimension(N+1), intent(in) :: theta
        real(kind=PR), dimension(N+1), intent(out) :: coord_x,coord_y

        coord_x=R*sin(theta*PI/180._PR) 
        coord_y=-R*cos(theta*PI/180._PR)

    end subroutine traject_cart 

    subroutine comparaison(to,T,N,thetai,thetai_c,err)
        integer, intent(in) :: N
        real(kind=PR), intent(in) :: to, T
        real(kind=PR), dimension(N+1), intent(in) :: thetai,thetai_c 
        real(kind=PR), dimension(N+1), intent(out) :: err 

        integer :: i 
        real(kind=PR) :: h 

        h=(T-to)/N

        do i=1,N+1
            err(i)=SQRT((thetai(i)-thetai_c(i))**2)
        end do 
    end subroutine comparaison 







    
    



end program 