MODULE LUSOLVER

  USE GLOBAL_DATA, ONLY : Ainv => Binv,maxsize

  IMPLICIT NONE
  REAL *8, DIMENSION (maxsize-1,maxsize-1)  :: u,l
  REAL *8, DIMENSION (maxsize-1)    :: b,y
!----------------------------------------------------------

CONTAINS

  SUBROUTINE LUINV(A,m,n)

    integer :: i,j,m,n
    real*8  :: A(m,n)
    real*8  :: den


    l = 1.0
    u = 0.0d0
    u(1,1) = a(1,1)
    u(1,2:n) = a(1,2:n)
    l(2:n,1) = a(2:n,1)/a(1,1)

    do i =2,n-1

       u(i,i) = a(i,i) - sum(l(i,1:i-1)*u(1:i-1,i))
       den = 1D0 / u(i,i)

       do j = i+1,n
          u(i,j) = a(i,j) - sum(l(i,1:i-1)*u(1:i-1,j))
          l(j,i) = (a(j,i) - sum(l(j,1:i-1)*u(1:i-1,i))) * den
       enddo
    enddo

    i=n
    u(i,i) = a(i,i) - sum(l(i,1:i-1)*u(1:i-1,i))


    do j = 1,n

       b(1:n) = 0.0d0
       b(j)   = 1.0d0

       y(1) = b(1)

       do i = 2,n
          y(i) = b(i) - sum(l(i,1:i-1)*y(1:i-1))
       enddo

       b(n) = y(n)/u(n,n)
       do i = n-1,1,-1
          b(i) = (y(i) - sum(u(i,i+1:n)*b(i+1:n)))/u(i,i)
       enddo

       Ainv(1:n,j) = b(1:n)

    enddo


    RETURN 
  END  SUBROUTINE LUINV

!********************************************************
SUBROUTINE LUSOLVE(A,f,m,n)

  integer :: i,j,m,n
  real*8  :: A(m,n),f(n)
  real*8  :: d


  l = 1.0
  u = 0.0d0
  u(1,1) = a(1,1)
  u(1,2:n) = a(1,2:n)
  l(2:n,1) = a(2:n,1)/a(1,1)

  do i =2,n-1

     u(i,i) = a(i,i) - sum(l(i,1:i-1)*u(1:i-1,i))
     d = 1.0 / u(i,i)

     do j = i+1,n
        u(i,j) = a(i,j) - sum(l(i,1:i-1)*u(1:i-1,j))
        l(j,i) = (a(j,i) - sum(l(j,1:i-1)*u(1:i-1,i))) * d
     enddo
  enddo

  i=n
  u(i,i) = a(i,i) - sum(l(i,1:i-1)*u(1:i-1,i))


  y(1) = f(1)

  do i = 2,n
     y(i) = f(i) - sum(l(i,1:i-1)*y(1:i-1))
  enddo

  f(n) = y(n)/u(n,n)
  do i = n-1,1,-1
     f(i) = (y(i) - sum(u(i,i+1:n)*f(i+1:n)))/u(i,i)
  enddo


  RETURN 
END SUBROUTINE LUSOLVE

END MODULE LUSOLVER
