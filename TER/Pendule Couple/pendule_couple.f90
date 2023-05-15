program couple 

    implicit none 

    integer, parameter :: PR=8
    real(PR), parameter :: PI=4._PR*ATAN(1._PR)
    integer :: N,i,meth
    real(kind=PR) :: m1,m2,g,l1,l2,C,to,T,thetao1,vto1,thetao2,vto2
    real(kind=PR), dimension(:), allocatable :: ti,theta1,pt1,theta2,pt2,coord_xt1,coord_yt1,coord_xt2,coord_yt2
    real(kind=PR), dimension(:), allocatable :: emi
    
   

    open(unit=1,file="parametres_3.dat",action="read")
    read(1,*) m1,m2,g,l1,l2,C,to,T,thetao1,vto1,thetao2,vto2,N,meth

    select case(meth)
    case(1)
        allocate(ti(N+1),theta1(N+1),pt1(N+1),theta2(N+1),pt2(N+1),coord_xt1(N+1),coord_yt1(N+1),coord_xt2(N+1),coord_yt2(N+1), &
        emi(N+1))
        call euler_implicite_nl(m1,m2,g,l1,l2,C,to,T,thetao1,vto1,thetao2,vto2,N,ti,theta1,pt1,theta2,pt2,emi)
        call traject_cart(l1,theta1,theta2,coord_xt1,coord_yt1,coord_xt2,coord_yt2)
        open(unit=3,file='euler_implicite1.dat',action='write')
        open(unit=4,file='euler_implicite2.dat',action='write')
        do i=1,N+1 
            write(3,*) ti(i),theta1(i)*180._PR/PI,pt1(i)/(m1*(l1**2)),coord_xt1(i),coord_yt1(i),emi(i)
            write(4,*) ti(i),theta2(i)*180._PR/PI,pt2(i)/(m2*(l2**2)),coord_xt2(i),coord_yt2(i),emi(i)
        end do
    case(2)
        allocate(ti(N+1),theta1(N+1),pt1(N+1),theta2(N+1),pt2(N+1),coord_xt1(N+1),coord_yt1(N+1),coord_xt2(N+1),coord_yt2(N+1), &
        emi(N+1))
        call euler_symplectique(m1,m2,g,l1,l2,C,to,T,thetao1,vto1,thetao2,vto2,N,ti,theta1,pt1,theta2,pt2,emi)
        call traject_cart(l1,theta1,theta2,coord_xt1,coord_yt1,coord_xt2,coord_yt2)
        open(unit=5,file='euler_symplectique1.dat',action='write')
        open(unit=6,file='euler_symplectique2.dat',action='write')
        do i=1,N+1 
            write(5,*) ti(i),theta1(i)*180._PR/PI,pt1(i)/(m1*(l1**2)),coord_xt1(i),coord_yt1(i),emi(i)
            write(6,*) ti(i),theta2(i)*180._PR/PI,pt2(i)/(m2*(l2**2)),coord_xt2(i),coord_yt2(i),emi(i)
        end do
    case(3)
        allocate(ti(N+1),theta1(N+1),pt1(N+1),theta2(N+1),pt2(N+1),coord_xt1(N+1),coord_yt1(N+1),coord_xt2(N+1),coord_yt2(N+1),&
        emi(N+1))
        call euler_explicite(m1,m2,g,l1,l2,C,to,T,thetao1,vto1,thetao2,vto2,N,ti,theta1,pt1,theta2,pt2,emi)
        call traject_cart(l1,theta1,theta2,coord_xt1,coord_yt1,coord_xt2,coord_yt2)
        open(unit=7,file='euler_explicite1.dat',action='write')
        open(unit=8,file='euler_explicite2.dat',action='write')
        do i=1,N+1 
            !write(7,*) ti(i),theta1(i)*180._PR/PI,pt1(i)/(m1*(l1**2)),coord_xt1(i),coord_yt1(i),emi(i)
            write(7,*) ti(i),emi(i)
            write(8,*) ti(i),theta2(i)*180._PR/PI,pt2(i)/(m2*(l2**2)),coord_xt2(i),coord_yt2(i),emi(i)
        end do
    case default 
        print*, "tapez 1"
    end select 


    
    
   
    
contains 

    subroutine LU(A,M)
        real(PR), dimension(:,:), allocatable, intent(in) :: A
        real(PR), dimension(:,:), allocatable, intent(out) :: M

        integer :: i,j,k

        M=A
        do k=1,size(A,1)-1
            do j=k+1,size(A,1)
                M(j,k)=M(j,k)/M(k,k)
            end do
            do j=k+1,size(A,1)
                do i=k+1,size(A,1)
                    M(j,i)=M(j,i)-M(j,k)*M(k,i)
                end do
            end do
        end do

    end subroutine

    subroutine reslu(M,b,x)
        real(PR), dimension(:,:), allocatable, intent(inout) :: M
        real(PR), dimension(:), allocatable, intent(in) :: b
        real(PR), dimension(:), allocatable, intent(inout) :: x

        real(PR), dimension(size(M,1),size(M,2)) :: L,U
        real(PR), dimension(size(M,1)) :: y
        integer :: i,j,k
        real(PR) :: S3,S4

        

        do k=1,size(M,1)
            L(k,k)=1
            U(k,k)=M(k,k)
            do j=1,k-1
                L(k,j)=M(k,j)
            end do
            do i=k+1,size(M,1)
                U(k,i)=M(k,i)
            end do
        end do

        S3=0

        y(1)=b(1)
        do i=2,size(M,1)
            do j=1,i-1
                S3=S3+L(i,j)*y(j)
            end do
            y(i)=b(i)-S3
            S3=0
        end do

        S4=0
        x(size(x))=y(size(x))/U(size(x),size(x))
        do i=size(M,1)-1,1,-1
            do j=i+1,size(M,1)
                S4=S4+U(i,j)*x(j)
            end do
            x(i)=(y(i)-S4)/U(i,i)
            S4=0
        end do

    end subroutine

    subroutine euler_implicite_nl(m1,m2,g,l1,l2,C,to,T,thetao1,vto1,thetao2,vto2,N,ti,theta1,pt1,theta2,pt2,emi)
        integer, intent(in) :: N !nombre de sous-intervalles 
        real(kind=PR), intent(in) :: m1,m2,g,l1,l2,C,to,T,thetao1,vto1,thetao2,vto2 
        real(kind=PR), dimension(N+1), intent(out) :: ti,theta1,pt1,theta2,pt2,emi

        integer :: i 
        real(kind=PR) :: h
        real(kind=PR), dimension(:), allocatable :: X,delta_X,FX 
        real(kind=PR), dimension(:,:), allocatable :: A,Id,J,V

        !X contient dans l'ordre theta1,theta2,p1 et p2
        !delta_X est mon inconnue : X(n+1)=X(n)+delta_X
        !FX=F(X) contient les 4 expressions des dérivées temporelles de theta1,theta2,pt1 et pt2
        !car on a des equations du type dtheta/dt=f(t,theta(t))
        allocate(X(4),delta_X(4),FX(4))
        !Id=identite, J est la jacobienne, A=Id-J et V me sert à faire LU 
        allocate(A(4,4),Id(4,4),J(4,4),V(4,4))

        !définition de mon pas de temps et initialisation 
        h=(T-to)/N 
        theta1(1)=thetao1
        theta2(1)=thetao2
        pt1(1)=m1*(l1**2)*vto1
        pt2(1)=m2*(l2**2)*vto2
        ti(1)=to
        emi(1)=((pt1(1)**2)/(m1*(l1**2))+(pt2(1)**2)/(m2*(l2**2)))/2._PR+g*(m1*l1*(1-COS(thetao1))+m2*l2*(1-COS(thetao2)))+ &
        (C*(thetao1-thetao2)**2)/2._PR

        !initialisation de Id,J et X
        Id=0._PR
        do i=1,4
            Id(i,i)=1._PR/h
        end do 
        
        J=0._PR
        J(1,3)=1/(m1*(l1**2))
        J(2,4)=1/(m2*(l2**2))
        J(3,2)=C
        J(4,1)=C

        X=(/theta1(1),theta2(1),pt1(1),pt2(1)/)

        do i=2,N+1
            !je change les valeurs de FX et certaines valeurs de J car elles dependent de theta1,2 et pt1,2
            FX(1)=pt1(i-1)/(m1*(l1**2))
            FX(2)=pt2(i-1)/(m2*(l2**2))
            FX(3)=-m1*g*l1*sin(theta1(i-1))-C*(theta1(i-1)-theta2(i-1))
            FX(4)=-m2*g*l2*sin(theta2(i-1))+C*(theta1(i-1)-theta2(i-1))
            J(3,1)=-m1*g*l1*cos(theta1(i-1))-C 
            J(4,2)=-m2*g*l2*cos(theta2(i-1))-C 

            ti(i)=ti(i-1)+h
            A=Id-J 
            !on résout le systeme A*delta_X=FX d'inconnue delta_X
            call LU(A,V)
            call reslu(V,FX,delta_X)
            !on calcule X en temps +1
            X=X+delta_X
            !je remplis mes differents vecteur pour pouvoir ensuite tracer ce dont j'ai envie
            theta1(i)=X(1)
            theta2(i)=X(2)
            pt1(i)=X(3)
            pt2(i)=X(4)
            !on calcule l'énergie avec les valeurs obtenues 
            emi(i)=((pt1(i)**2)/(m1*(l1**2))+(pt2(i)**2)/(m2*(l2**2)))/2._PR+g*(m1*l1*(1-COS(theta1(i)))+m2*l2*(1-COS(theta2(i))))&
            +(C*(theta1(i)-theta2(i))**2)/2._PR
        end do 

    end subroutine euler_implicite_nl 

    subroutine euler_symplectique(m1,m2,g,l1,l2,C,to,T,thetao1,vto1,thetao2,vto2,N,ti,theta1,pt1,theta2,pt2,emi)
        integer, intent(in) :: N !nombre de sous-intervalles 
        real(kind=PR), intent(in) :: m1,m2,g,l1,l2,C,to,T,thetao1,vto1,thetao2,vto2 
        real(kind=PR), dimension(N+1), intent(out) :: ti,theta1,pt1,theta2,pt2,emi

        integer :: i 
        real(kind=PR) :: h

        !définition de mon pas de temps et initialisation 
        h=(T-to)/N 
        theta1(1)=thetao1
        theta2(1)=thetao2
        pt1(1)=m1*(l1**2)*vto1
        pt2(1)=m2*(l2**2)*vto2
        ti(1)=to
        emi(1)=((pt1(1)**2)/(m1*(l1**2))+(pt2(1)**2)/(m2*(l2**2)))/2._PR+g*(m1*l1*(1-COS(thetao1))+m2*l2*(1-COS(thetao2)))+ &
        (C*(thetao1-thetao2)**2)/2._PR


        do i=2,N+1
            ti(i)=ti(i-1)+h
            theta1(i)=theta1(i-1)+h*pt1(i-1)/(m1*(l1**2))
            theta2(i)=theta2(i-1)+h*pt2(i-1)/(m2*(l2**2))
            pt1(i)=pt1(i-1)-h*(g*m1*l1*SIN(theta1(i))+C*(theta1(i)-theta2(i)))
            pt2(i)=pt2(i-1)-h*(g*m2*l2*SIN(theta2(i))-C*(theta1(i)-theta2(i)))
            emi(i)=((pt1(i)**2)/(m1*(l1**2))+(pt2(i)**2)/(m2*(l2**2)))/2._PR+g*(m1*l1*(1-COS(theta1(i)))+m2*l2*(1-COS(theta2(i))))&
            +(C*(theta1(i)-theta2(i))**2)/2._PR
        end do 



    end subroutine euler_symplectique 

    subroutine euler_explicite(m1,m2,g,l1,l2,C,to,T,thetao1,vto1,thetao2,vto2,N,ti,theta1,pt1,theta2,pt2,emi)
        integer, intent(in) :: N !nombre de sous-intervalles 
        real(kind=PR), intent(in) :: m1,m2,g,l1,l2,C,to,T,thetao1,vto1,thetao2,vto2 
        real(kind=PR), dimension(N+1), intent(out) :: ti,theta1,pt1,theta2,pt2,emi

        integer :: i 
        real(kind=PR) :: h

        !définition de mon pas de temps et initialisation 
        h=(T-to)/N 
        theta1(1)=thetao1
        theta2(1)=thetao2
        pt1(1)=m1*(l1**2)*vto1
        pt2(1)=m2*(l2**2)*vto2
        ti(1)=to
        emi(1)=((pt1(1)**2)/(m1*(l1**2))+(pt2(1)**2)/(m2*(l2**2)))/2._PR+g*(m1*l1*(1-COS(thetao1))+m2*l2*(1-COS(thetao2)))+ &
        (C*(thetao1-thetao2)**2)/2._PR

        do i=2,N+1
            ti(i)=ti(i-1)+h
            theta1(i)=theta1(i-1)+h*pt1(i-1)/(m1*(l1**2))
            theta2(i)=theta2(i-1)+h*pt2(i-1)/(m2*(l2**2))
            pt1(i)=pt1(i-1)-h*(g*m1*l1*SIN(theta1(i-1)*PI/180._PR)+C*(theta1(i-1)-theta2(i-1)))
            pt2(i)=pt2(i-1)-h*(g*m2*l2*SIN(theta2(i-1)*PI/180._PR)-C*(theta1(i-1)-theta2(i-1)))
            emi(i)=((pt1(i)**2)/(m1*(l1**2))+(pt2(i)**2)/(m2*(l2**2)))/2._PR+g*(m1*l1*(1-COS(theta1(i)))+m2*l2*(1-COS(theta2(i))))&
            +(C*(theta1(i)-theta2(i))**2)/2._PR
        end do 



    end subroutine euler_explicite

    subroutine traject_cart(R,theta1,theta2,coord_xt1,coord_yt1,coord_xt2,coord_yt2)
        ! Trace la trajectoire de l'objet dans un repère cartésien (initialement polaire).
        real(kind=PR), intent(in) :: R
        real(kind=PR), dimension(N+1), intent(in) :: theta1,theta2
        real(kind=PR), dimension(N+1), intent(out) :: coord_xt1,coord_yt1,coord_xt2,coord_yt2

        integer :: i 

        do i=1,N+1
            coord_xt1(i)=R*sin(theta1(i)*PI/180._PR) 
            coord_yt1(i)=-R*cos(theta1(i)*PI/180._PR)
            coord_xt2(i)=R*sin(theta2(i)*PI/180._PR) 
            coord_yt2(i)=-R*cos(theta2(i)*PI/180._PR)
        end do 
    end subroutine traject_cart 

end program 