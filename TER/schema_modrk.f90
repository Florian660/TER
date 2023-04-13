module schema_mod
    use constantes_mod
    use edo_mod
    implicit none

    type :: rk_type
        integer                               :: s
        real(PR), dimension(:,:), allocatable :: A
        real(PR), dimension(:), allocatable   :: b
    end type

contains

    subroutine un_pas_temps(u,edo,rk, dt)
        real(PR), dimension(:), intent(inout) :: u
        type(edo_type), intent(inout)         :: edo
        type(rk_type), intent(inout)          :: rk
        real(PR), intent(in)                  :: dt
        integer                               :: i,j
        real(PR),dimension(:,:),allocatable   :: k
        real(PR),dimension(:),allocatable     :: v,som1, som2


        allocate(k(size(u),rk%s))
        allocate(v(size(u)),som1(size(u)))
        allocate(som2(size(u)))

        k(:,1) = fh(u)
        if (rk%s >1)then
            do i = 2, rk%s
                som1 = 0._PR
                do j = 1,i-1
                    som1 = som1 + rk%A(i,j)*k(:,j)
                end do
                v = u + dt*som1
                k(:,i) = fh(v)
            end do
        end if

        do i = 1, size(u)
            som2 = 0._PR
            do j = 1,rk%s
                som2 = som2 + rk%b(j)*k(:,j)
            end do
            u(i) = u(i)+ dt*som2(i)
        end do

        deallocate(v,k)
        deallocate(som1,som2)
    end subroutine

    subroutine init_rk(rk,edo)
        type(rk_type), intent(inout) :: rk
        type(edo_type), intent(in)   :: edo

        rk%s = edo%sch
        Selectcase (edo%sch)
            case (1)
                allocate(rk%A(1,0),rk%b(1))
                rk%A = 0._PR
                rk%b = 1._PR

            case (2)
                allocate(rk%A(2,2),rk%b(2))
                rk%A = 0._PR
                rk%A(2,1) = 1._PR
                rk%b = 1/2._PR

            case(3)
                allocate(rk%A(3,3),rk%b(3))
                rk%A = 0._PR
                rk%A(2,1) = 1/2._PR
                rk%A(3,1) = -1._PR
                rk%A(3,2) = 2._PR
                rk%b = 1/6._PR
                rk%b(2) = 2/3._PR

            case(4)
                allocate(rk%A(4,4),rk%b(4))
                rk%A = 0._PR
                rk%A(2,1) = 1/2._PR
                rk%A(3,2) = 1/2._PR
                rk%A(4,3) = 1._PR
                rk%b = 1/6._PR
                rk%b(2) = 1/3._PR
                rk%b(3) = 1/3._PR
        end select

    end subroutine

    subroutine free_rk(rk)
        type(rk_type), intent(inout) :: rk

        deallocate(rk%A,rk%b)
    end subroutine





    !
    ! subroutine calcul_edo1(u,edo,rk)
    !     type(edo_type), intent(inout)         :: edo
    !     type(rk_type),intent(inout)           :: rk
    !     real(PR), dimension(:), intent(inout) :: u
    !     real(PR)                              :: t
    !
    !     t = 0._PR
    !
    !     do while(t<edo%tfin)
    !         call un_pas_temps(u,edo,rk)
    !
    !         t = t+edo%deltat
    !
    !     end do
    !     print*, u-val_exa
    ! end subroutine
end module schema_mod
