module constantes_mod
    implicit none
    integer,parameter   :: PR=8
    real(PR), parameter :: val_exa = 2._PR !valeur exacte pour edo1 a t = 1
    real(PR), parameter :: g=9.81, l = 0.9, w= sqrt(g/l), m = 0.5

contains
    subroutine printMat(M)
            real(PR), dimension(:,:), intent(in)  :: M
            integer                               :: i

            do i = 1,size(M,1)
                print*,M(i,:)

            end do
        end subroutine

        function threel(th0,t)result(res)
            real(PR), intent(in) :: th0,t
            real(PR)             :: res


            res = th0*cos(w*t)

        end function threel

        function vitreel (th0,t) result(res)
            real(PR), intent(in) :: th0,t
            real(PR)             :: res

            res = - th0*w*sin(w*t)


        end function vitreel
end module constantes_mod
