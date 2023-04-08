program pendules

    implicit none 

    integer, parameter :: PR=8
    real(PR), parameter :: PI=4*ATAN(1._PR)
    integer :: N,i,j
    real(kind=PR) :: to,thetao,T,g,l,m,vo,h,S
    real(kind=PR), dimension(:), allocatable :: ti,thetai,vi,coord_x, coord_y, thetai_c, vi_c,coord_x_c,coord_y_c, err
    real(kind=PR), dimension(:), allocatable :: emi,emic

   

 

    open(unit=1,file="parametres.dat",action="read")
    read(1,*) m,g,l,to,T,thetao, vo, N

    !l'argument dans un cosinus ou sinus est de base en radian

    
    ! CALCUL DE L'ERREUR 
    allocate(err(8))
    
    open(unit=3,file='erreur.dat',action="write")
    do i=1,8  ! on fait varier le pas de temps de 10^-3 à 10^-5 
        S=0
        N=50*2**(i+4) !pas de temps de 2.5*10**(-2) à 2.4*10**(-5)
       
        
        allocate(ti(N+1),thetai(N+1),vi(N+1),coord_x(N+1),coord_y(N+1), thetai_c(N+1),vi_c(N+1),coord_x_c(N+1), &
        coord_y_c(N+1),emi(N+1), emic(N+1))  
        call realite(g,l,to,thetao,vo,T,N,ti,thetai_c,vi_c,m,emic)  !on calcule la "vraie" valeur de theta pour chaque pas de temps
        
        call euler_explicite(m,g,l,to,thetao,vo,T,N,ti,thetai,vi,emi) ! de même pour le schéma d'euler 
       
   
        do j=1,N
            S=S+(thetai_c(j)-thetai(j))**2
        end do 
        h=(T-to)/N
        err(i)=SQRT(S*h)  ! on calcule l'erreur en norme L2
           
       
        

        write(3,*) h, err(i)
        deallocate(ti,thetai,vi,coord_x,coord_y, thetai_c, vi_c, coord_x_c, coord_y_c,emi,emic)
    end do 
    
    deallocate(err)

    
    ! allocate(ti(N+1),thetai(N+1),vi(N+1),coord_x(N+1),coord_y(N+1), thetai_c(N+1),vi_c(N+1),coord_x_c(N+1), &
    !     coord_y_c(N+1),err(N+1),emi(N+1),emic(N+1))
    ! call euler_explicite(m,g,l,to,thetao,vo,T,N,ti,thetai,vi,emi)
    ! call traject_cart(l,thetai,coord_x,coord_y)
    ! call realite(g,l,to,thetao,vo,T,N,ti,thetai_c,vi_c,m,emic)
    ! call traject_cart(l,thetai_c,coord_x_c,coord_y_c)
    ! call comparaison(to,T,N,thetai,thetai_c,err)

    ! open(unit=2,file="euler.dat",action="write")

    ! do i=1,N+1
    !     write(2,*) ti(i), thetai(i), vi(i), coord_x(i), coord_y(i), thetai_c(i), vi_c(i), coord_x_c(i), &
    !     coord_y_c(i), err(i)
      
    ! end do 
    ! deallocate(ti,thetai,vi,coord_x,coord_y, thetai_c, vi_c, coord_x_c, coord_y_c,err,emi,emic)
   

    ! close(1)
    ! close(2)

    !ENERGIE MECANIQUE 

    ! allocate(ti(N+1),thetai(N+1),vi(N+1),thetai_c(N+1),vi_c(N+1),emi(N+1),emic(N+1))
    ! call euler_explicite(m,g,l,to,thetao,vo,T,N,ti,thetai,vi,emi)
    ! call realite(g,l,to,thetao,vo,T,N,ti,thetai_c,vi_c,m,emic)
    ! open(10,file='energie.dat',action="write")
  
    ! do i=1,N+1
    !     write(10,*) ti(i),emi(i),emic(i)
        
    ! end do 

    ! deallocate(ti,thetai,vi,thetai_c,vi_c,emi,emic)
    ! close(10)
    



    contains

    function f (a,b,c) result(y)
        real(kind=PR), intent(in) :: a,b,c

        real(kind=PR) :: y 

        y=-(a/b)*c

    end function f 

    function em(a,b,c,d,e) result(y)
        real(kind=PR), intent(in) :: a,b,c,d,e

        real(kind=PR) :: y 

        y=(1/2._PR)*a*(b**2)*(c**2)-a*b*d*(1-(e**2)/2._PR)

    end function em


   
    subroutine euler_explicite (m,g,l,to,thetao,vo,T,N,ti,thetai,vi,emi)
        integer, intent(in) :: N !nombre de sous-intervalles 
        real(kind=PR), intent(in) :: m,g, l, to, thetao, T, vo 
        
        real(kind=PR), dimension(N+1), intent(out) :: ti,thetai,vi,emi

        integer :: i 
        real(kind=PR) :: h

        

        h=(T-to)/N 
        thetai(1)=thetao
        vi(1)=vo
        ti(1)=to
        emi(1)=em(m,l,vo,g,thetao)

        do i=2,N+1
            ti(i)=ti(i-1)+h
            thetai(i)=thetai(i-1)+h*vi(i-1)
            vi(i)=vi(i-1)+h*f(g,l,thetai(i-1))
            emi(i)=em(m,l,vi(i),g,thetai(i))
        end do 

    end subroutine euler_explicite 

    subroutine euler_implicite(g,l,to,thetao,vo,T,N,ti,thetai,vi)   !petites oscillations
        integer, intent(in) :: N !nombre de sous-intervalles 
        real(kind=PR), intent(in) :: g, l, to, thetao, T, vo 
        
        real(kind=PR), dimension(N+1), intent(out) :: ti,thetai,vi

        integer :: i 
        real(kind=PR) :: h,p,r

        p=l+g*(h**2)
        r=g*(h**2)

        h=(T-to)/N 
        thetai(1)=thetao
        vi(1)=vo
        ti(1)=to

        
        do i=2,N+1
            ti(i)=ti(i-1)+h
            thetai(i)=(l/p)*(thetai(i-1)+h*vi(i-1))
            vi(i)=vi(i-1)*(1-(r/p))-((h*g)/p)*thetai(i-1)
        end do 
    
    end subroutine euler_implicite

    subroutine realite(g,l,to,thetao,vo,T,N,ti,thetai,vi,m,emic)
        integer, intent(in) :: N !nombre de sous-intervalles 
        real(kind=PR), intent(in) :: g, l, to, thetao, T, vo,m 
        real(kind=PR), dimension(N+1), intent(out) :: ti,thetai,vi,emic

        integer :: i 
        real(kind=PR) :: h,omega

        h=(T-to)/N 
        thetai(1)=thetao
        ti(1)=to
        vi(1)=vo
        emic(1)=em(m,l,vo,g,thetao)
       

        omega=SQRT(g/l)

        do i=2,N+1
            ti(i)=ti(i-1)+h
            thetai(i)=thetao*COS(omega*ti(i))+(vo/omega)*SIN(omega*ti(i))
            vi(i)=-thetao*omega*SIN(omega*ti(i))+vo*COS(omega*ti(i))
            emic(i)=em(m,l,vi(i),g,thetai(i))
            
        
        end do

    end subroutine realite


    subroutine traject_cart(R, theta, coord_x, coord_y)
        ! Trace la trajectoire de l'objet dans un repère cartésien (initialement polaire).
        real(kind=PR), intent(in) :: R
        real(kind=PR), dimension(N+1), intent(in) :: theta
        real(kind=PR), dimension(N+1), intent(out) :: coord_x, coord_y

        
        coord_x = R*sin(theta*PI/180._PR) 
        coord_y = -R*cos(theta*PI/180._PR)

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