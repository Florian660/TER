!programme principal
program eulerimplicit
  use fonctions
  implicit none  
  real*8,dimension(:,:),allocatable::A,L,U,mat
  real*8,dimension(:),allocatable::Un,Un1,y
  real*8::Beta1,alpha,dt,dx,tf,D
  integer::i,j,nx,nt
  D=0.01
  tf=1.0D0
  dt=0.1d0
  dx=0.1d0
 
  nx=floor(1./dx)
  nt=floor(tf/dt)
 
  allocate(A(nx,nx))
  allocate(L(nx,nx))
  allocate(U(nx,nx))
  allocate(y(nx))
  allocate(Un(nx))
  allocate(Un1(nx))
  allocate(mat(nx,nt))
  Beta1=(1+2*(dt/dx**2))
  alpha=dt/dx**2
  A(1,1)=1.d0
  A(nx,nx)=1.d0
  Un(1)=0.d0
  Un(nx)=0.d0
 
  do i=2,nx-1
     A(i,i)=Beta1
     A(i,i+1)=-alpha
     A(i,i-1)=-alpha
  end do
 
 
  Un(1)=0.d0
  Un(nx)=0.d0
  do j=2,nx-1
     Un(j)=5.d0
  end do
  open(unit=24,file='resultatsimplicistp')
  print*,A
  do i=1,nt
     call LU(nt,D*A,L,U) !on multiplie la matrice A par le coefficient D'
     call Resolution(nt,L,U,Un,Un1)
     mat(:,i) = Un1
     Un=Un1
     ! print*,'la matrice L est :',L
     ! print*,'la matrice U est :',U
     print*, 'la matrice Un1 est :',Un1
  end do
  call writegnuplot('implicit.dat',mat)
  close(24)
 
  deallocate(A)
  deallocate(L)
  deallocate(U)
  deallocate(y)
  deallocate(Un)
  deallocate(Un1)
 
end program eulerimplicit
