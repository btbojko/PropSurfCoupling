! ********************************************************************
! ********************************************************************
! Updated: TLJ; 12/19/2016
! Filename: pyrolysis.f90
! ********************************************************************
! ********************************************************************
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
!
!     Subroutine computes the local burn rate from
!       a pyrolysis law
!
!     subroutines: PYROLYSIS
!
! ********************************************************************
SUBROUTINE PYROLYSIS

  USE GLOBAL_DATA

  IMPLICIT NONE

! Local Variables
  INTEGER :: i, k
  REAL*8  :: phixAVG,phizAVG,rbxavg,rbzavg
!--------------------------------------------------------------------

  DO k=drange(3)-1,drange(4)+1
  DO i=drange(1)-1,drange(2)+1

     rb(i,k) = da_s(i,k,0)*exp(-theta_s(i,k,0)/f(i,k,0,1))
     d_rb(i,k) = theta_s(i,k,0)/(f(i,k,0,1)**2)*rb(i,k)

     ft0(i,k) = rb(i,k)*sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)
     d_ft0(i,k) = d_rb(i,k)*sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)

     vvel(i,k,0) = rhos(i,k,0)*ft0(i,k)
     d_vvel(i,k,0) = rhos(i,k,0)*d_ft0(i,k)

     rbxavg = half*(rb(i,k)+rb(i+1,k))
     rbzavg = half*(rb(i,k)+rb(i,k+1))
     phizAVG = half*(dphidz(i,k)+dphidz(i+1,k))
     phixAVG = half*(dphidx(i,k)+dphidx(i,k+1))

     vec(i,k,0) = rb(i,k)*sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)
     vec(i,k,1) = rbxavg*sqrt(1.0d0+ dphidxa(i,k)**2 + phizAVG**2)
     vec(i,k,2) = rbzavg*sqrt(1.0d0+ phixAVG**2 + dphidza(i,k)**2)

  END DO
  END DO
  
  RETURN
END SUBROUTINE PYROLYSIS
!*********************************************************************
