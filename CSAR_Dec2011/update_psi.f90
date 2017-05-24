! Heterogeneous Propellant 3D Combustion Code
! Constant Mass Flux
! DOE-ASCI Center for Simulation of Advanced Rockets
! University of Illinois at Urbana-Champaign
! Author: Thomas L. Jackson
! Date last Modified: 16 January 2001
! Parallelized by M. Campbell on May 9, 2001
! //////////////////////////////////////////////////
!
! ********************************************************************
! PROPRIETARY CODE - DO NOT USE WITHOUT PERMISSION OF CSAR!
! ********************************************************************
!
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
!
!     Level Set
!
!     Subroutine updates the level-set function psi(i,j)
!       used in the solid to distinguish between pure
!       AP and binder; also updates surface and solid
!       values
!
!     Solid level-set function to compute solid values
!        psi = +1    AP, oxidizer
!        psi = -1    Binder, fuel
!
! ********************************************************************
SUBROUTINE UPDATE_PSI(t,tout)

USE GLOBAL_DATA

IMPLICIT NONE

INTEGER ::  i ,n, j, k, m, kkk
INTEGER :: yrange,zmin,zmax,xmin,xmax
REAL*8  :: xcut,ycut,t,tout,ycutold
!
! computes the level set for each mode separately
!
k = 0
if (ipack == 1) THEN
    DO j=drange(5),drange(6)
        DO i=drange(1),drange(2) 
            xcut = x(i)
            ycut = phi(i,k) - y(j)
            ycut = mod(ycut,period)
            call levelset_2D(xcut,ycut,i,j)    
        END DO
    END DO
else if (ipack == 0) Then
    do i = drange(1)-1,drange(2)+1
        psi(i,:,:) = one          !AP
        if (abs(x(i)) > xloc) psi(i,:,:) = -one      !BINDER
    enddo
else
     do i = drange(1)-1,drange(2)+1
        psi(i,:,:) = one          !AP
        if (alpha_pack  /= 1) psi(i,:,:) = -one      !BINDER
     enddo   
endif

CALL FILL_PSI_GHOST
!  if(ipack >= 0) then
CALL SMOOTH
!  else
CALL SOLID_PROP
!  endif

CALL PYROLYSIS

!  CALL VISUALIZE_FIELD 
  
!  CALL PRN_TMP
  
!------------------------------------------------------------------------
RETURN
END SUBROUTINE UPDATE_PSI
