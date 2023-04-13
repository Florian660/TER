module edo_mod
    use constantes_mod
    implicit none
    type :: edo_type
        integer                             :: dim
        real(PR)                            :: tfin
        real(PR), dimension(:), allocatable :: u0
        integer                             :: sch, n
    end type

contains


    function f(v) result(res)  ! cas vect a deux composantes
        real(PR), dimension(2), intent(in) :: v
        real(PR), dimension(2)            :: res


        res(1) = -(g/l)*v(2)
        res(2) = v(1)


    end function f

    function fh(v) result(res)
        real(PR), dimension (2), intent(in) :: v
        real(PR), dimension (2)             :: res
        ! v = (p theta)
        res(1) = -(g/l)*SIN(v(2))
        res(2) = v(1)

    end function fh



    subroutine init_edo(a)
        type(edo_type), intent(inout) :: a

        open(unit=1, file='parametresrk.dat', form='formatted')

        read(1,*) a%dim
        read(1,*) a%n
        read(1,*) a%tfin
        read(1,*) a%sch

        allocate(a%u0(a%dim))

        read(1,*) a%u0 !u = (v,th)

        close(1)


    end subroutine


    subroutine free_edo(a)
        type(edo_type), intent(inout) :: a

        deallocate(a%u0)
    end subroutine

end module edo_mod
