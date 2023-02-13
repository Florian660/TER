program hello  

    implicit none 

    integer, parameter :: PR=8
    integer :: N,i
    real(kind=PR) :: to,thetao,T,g,l,vo
    real(kind=PR), dimension(:), allocatable :: ti,thetai,vi

   

    

    open(unit=1,file="parametres.dat",action="read")
    read(1,*) g,l,to,T,thetao, vo, N

    allocate(ti(N+1),thetai(N+1),vi(N+1))

    

    

    call euler(g,l,to,thetao,vo,T,N,ti,thetai,vi)

    open(unit=2,file="euler.dat",action="write")

    do i=1,N+1
        write(2,*) ti(i), thetai(i), vi(i)
    end do 

   

   

    

    deallocate(ti,thetai,vi)

    close(1)
    close(2)

    contains

    function f (a,b,c) result(y)
        real(kind=PR), intent(in) :: a,b,c

        real(kind=PR) :: y 

        y=-(a/b)*c

    end function f 

    subroutine euler (g,l,to,thetao,vo,T,N,ti,thetai,vi)
        integer, intent(in) :: N !nombre de sous-intervalles 
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

    
    



end program 