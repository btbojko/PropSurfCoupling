!**************************************************************************
SUBROUTINE LEVELSET_2D(xcut,ycut,i,j)

  USE GLOBAL_DATA

  IMPLICIT NONE

!  Local Variables
  INTEGER :: n,i,j,k 
  REAL*8  :: psinew, r2, psi1
  REAL*8  :: xcut, ycut, a_1,a_2, theta,a_3
!-----------------------------------------------------------------

! computes the level set for each mode separately

  k = 0
  psinew = 100.0d0
  do n=1,ncircles
!-modi0
     a_1 = ycut - ycenter(n)
     a_2 = xcut - xcenter(n)
     a_3 = a_1/a_2
     if(a_2 == 0.0d0)then
        theta = 0.0d0
     else
        theta = ATAN(a_3)
     endif
     r2 = (alpha_test*rad(n)*cos(theta))**2+(1/alpha_test*rad(n)*sin(theta))**2   


!     r2 = rad(n)*rad(n)
     psi1 = ((xcut-xcenter(n))**2 &
          + (ycut-ycenter(n))**2) - r2
     psinew = min(psi1,psinew)
  end do
  psi(i,k,j) = -psinew*1.d4

!-modi0
  if (psi(i,k,j) .gt. 0.d0) psi(i,k,j) = .25
!  if (psi(i,k,j) .lt. 0.d0) psi(i,k,j) = -1.

  RETURN
END SUBROUTINE LEVELSET_2D
!******************************************************************
