program double 

    implicit none 

    integer, parameter :: PR=8
    real(PR), parameter :: PI=4._PR*ATAN(1._PR)
    integer :: N,i
    real(kind=PR) :: to,thetao,T,g,l,m,vto,psio,vpo
    real(kind=PR), dimension(:), allocatable :: ti,thetai,vti,psi,vpi,coord_xt,coord_yt,coord_xp,coord_yp
   

    open(unit=1,file="parametres_2.dat",action="read")
    read(1,*) m,g,l,to,T,thetao,vto,psio,vpo,N

    allocate(ti(N+1),thetai(N+1),vti(N+1),psi(N+1),vpi(N+1),coord_xt(N+1),coord_yt(N+1),coord_xp(N+1),coord_yp(N+1))
    call euler_explicite(m,g,l,to,T,thetao,vto,psio,vpo,N,ti,thetai,vti,psi,vpi)
    call traject_cart(l,thetai,psi,coord_xt,coord_yt,coord_xp,coord_yp)
    open(unit=4,file='euler_explicite.dat',action='write')
    
    do i=1,N+1 
        write(4,*) ti(i),thetai(i)*180._PR/PI,vti(i),psi(i)*180._PR/PI,vpi(i), &
        coord_xt(i),coord_yt(i),coord_xp(i),coord_yp(i)
    end do

    



contains 
    subroutine euler_explicite (m,g,l,to,T,thetao,vto,psio,vpo,N,ti,thetai,vti,psi,vpi) 
        integer, intent(in) :: N !nombre de sous-intervalles 
        real(kind=PR), intent(in) :: m,g,l,to,T,thetao,vto,psio,vpo
        real(kind=PR), dimension(:), allocatable, intent(inout) :: ti,thetai,vti,psi,vpi

        integer :: i 
        real(kind=PR) :: h

        h=(T-to)/N 
        ti(1)=to
        thetai(1)=thetao
        vti(1)=vto
        psi(1)=psio
        vpi(1)=vpo
        !m1=m2 l1=l2

        do i=1,N
            ti(i+1)=ti(i)+h
            thetai(i+1)=thetai(i)+h*vti(i)
            thetai(i+1)=MODULO(thetai(i+1),2._PR*PI)
            vti(i+1)=vti(i)+h*(-2._PR*m*g*SIN(thetai(i))+m*(g*COS((psi(i)-thetai(i)))*SIN(psi(i)) &
            +l*(vti(i)**2)*SIN(psi(i)-thetai(i))*COS(psi(i)-thetai(i))+ &
            l*(vpi(i)**2)*SIN(psi(i)-thetai(i))))/(2*m*l-m*l*(COS(psi(i)-thetai(i)))**2)
            psi(i+1)=psi(i)+h*vpi(i)
            psi(i+1)=MODULO(psi(i+1),2._PR*PI)
            vpi(i+1)=vpi(i)+h*(-2*m*(g*SIN(thetai(i))*COS(psi(i)-thetai(i))-g*SIN(psi(i))- &
            l*(vti(i)**2)*SIN(psi(i)-thetai(i)))+m*l*(vpi(i)**2)*COS(psi(i)-thetai(i))* &
            SIN(psi(i)-thetai(i)))/(m*l*(COS(psi(i)-thetai(i)))**2-2*m*l)
        end do 

       
    end subroutine euler_explicite

    subroutine traject_cart(R,theta,psi,coord_xt,coord_yt,coord_xp,coord_yp)
        ! Trace la trajectoire de l'objet dans un repère cartésien (initialement polaire).
        real(kind=PR), intent(in) :: R
        real(kind=PR), dimension(N+1), intent(in) :: theta,psi
        real(kind=PR), dimension(N+1), intent(out) :: coord_xt,coord_yt,coord_xp,coord_yp

        integer :: i 

        do i=1,N+1
            coord_xt(i)=R*sin(theta(i)) 
            coord_yt(i)=-R*cos(theta(i))
            coord_xp(i)=coord_xt(i)+R*sin(psi(i)) 
            coord_yp(i)=coord_yt(i)-R*cos(psi(i))
        end do 

    end subroutine traject_cart 

end program 