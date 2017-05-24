
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
!
!     Subroutine computes the local burn rate from
!       a pyrolysis law

! ********************************************************************
SUBROUTINE PYROLYSIS

  USE GLOBAL_DATA

  IMPLICIT NONE

! Local Variables
  INTEGER ::  i, k
!---------------------------------------------------------

    DO k=drange(3)-1,drange(4)+1
       DO i=drange(1)-1,drange(2)+1
          rb(i,k) = da_s(i,k,0)*exp(-theta_s(i,k,0)/f(i,k,0,1))
          d_rb(i,k) = theta_s(i,k,0)/f(i,k,0,1)**2*rb(i,k)
          ft0(i,k) = rb(i,k)*sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)      !evaluated twice
          d_ft0(i,k) = d_rb(i,k)*sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)
       END DO
    END DO

  
 
  
  RETURN
END SUBROUTINE PYROLYSIS
!***********************************************************
