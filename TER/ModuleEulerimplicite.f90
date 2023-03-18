module fonctions
 
  !module principal contenant les subroutines
  implicit none
 
contains
 
!**************************************************************************
!Résolution d'un système linéaire triangulaire supérieur
 
subroutine Resol_Triang_Sup(n,A,B,X)
 
  implicit none
 
  integer, intent(in) :: n
  real*8, dimension(n,n),intent(in) :: A
  real*8, dimension(n), intent(in) :: B
  real*8, dimension(n), intent(out) :: X
  integer :: i, j
  real*8 :: somme
 
 
  !initialisation de X et somme
  X = 0
  somme = 0
 
  !algorithme de résolution
  X(n) = B(n)/A(n,n)
 
  Do i = n-1,1,-1
     somme = B(i)
     Do j = n, i, -1 
        somme = somme - A(i,j)*X(j)
     End Do
     X(i) = somme/A(i,i)
  End Do
 
end subroutine Resol_Triang_Sup
 
!*****************************************************************************
!Résolution d'un système linéaire triangulaire inférieur
 
subroutine Resol_Triang_Inf(n,A,B,X)
 
  implicit none
 
  integer, intent(in) :: n
  real*8, dimension(n,n),intent(in) :: A
  real*8, dimension(n), intent(in) :: B
  real*8, dimension(n), intent(out) :: X
  integer :: i, j
  real*8 :: somme
 
 
  !initialisation de X et somme
  X = 0
  somme = 0
 
  !algorithme de résolution
  X(1) = B(1)/A(1,1)
 
  Do i = 2, n
     somme = B(i)
     Do j = 1, i-1
        somme = somme - A(i,j)*X(j)
     End Do
     X(i) = somme/A(i,i)
  End Do
 
end subroutine Resol_Triang_Inf
!**************************************************************************
 
!Factorisation LU
 
subroutine LU(n,A,L,U)
 
  implicit none
 
  integer, intent(in) :: n
  real*8, dimension(n,n), intent(in) :: A
  real*8, dimension(n,n), intent(out) :: L, U
  integer :: i, j, k
 
  U = A
  L = 0
 
  Do i = 1, n
     L(i,i) = 1
  End Do
 
  Do i = 1, n
     Do j = i+1, n
        L(j,i) = U(j,i)/U(i,i)
        U(j,i) = 0
        Do k = i+1, n
           U(j,k) = U(j,k) - L(j,i)*U(i,k)
        End Do
     End Do
  End Do
 
end subroutine LU
!***************************************************************************
!Resolution LU
subroutine Resolution(n,L,U,B,X)
 
  implicit none
 
  integer, intent(in) :: n
  real*8, dimension(n,n), intent(in) :: L, U
  real*8, dimension(n), intent(in) :: B
  real*8, dimension(n), intent(out) :: X
  real*8, dimension(n) :: Y
  integer :: i, j
 
  call Resol_Triang_Inf(n,L,B,Y)
  call Resol_Triang_Sup(n,U,Y,X)
 
 
end subroutine Resolution
 
!***************************************************************
 
!ecriture dans un fichier
 
subroutine WriteGnuPlot(nom_fichier, v)
    implicit none
    character(len=*) :: nom_fichier
    real*8, dimension(:,:) :: v
    integer :: i,j
    open(unit = 11, file = nom_fichier, form = 'formatted', &
         status = 'unknown', action = 'write')
    do i = 1, size(v,1)
       do j=1,size(v,2)
          write(11, * )i, j, v(i,j)
       end do
       write(11, *)
    end do
    close(11)
  end subroutine WriteGnuPlot
 
end module fonctions
