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
!     Grid Generation
!
!     This function establishes the computational grid.
!
! ********************************************************************
SUBROUTINE GRIDSETUP
 
USE data_types
USE GLOBAL_DATA
  
IMPLICIT NONE

!-------------------------------------------------------------------------
!     Local Variables

TYPE(mesh_pointers), POINTER :: current_mesh 
INTEGER ::  i, j, k , jj, mny, ye
REAL*8  :: den, der1, der2
REAL*8  :: eta, eta2
REAL*8  :: dyf, dyc
!-------------------------------------------------------------------------

DO i = drange(1)-1, drange(2) +1
    x(i) = xstart_MG + dfloat(i+ddrange(1))*dx
END DO
  
! z-grid; uniform and with no mapping
DO k = drange(3)-1, drange(4)+1
    z(k) = zero
END DO
  
! y-grid  2 loops becasue of staggering
DO j = 0, ny
    eta  = dfloat(j+ddrange(5))*dy
    eta2 = (eta/yend)**2
    den = 2.0d0 - eta2
    y(j) = eta/den**c1y
    der1 = 1.0d0/(den**c1y)+2.0d0*eta2*c1y/den**(c1y+1.0d0)
    der2 = 2.0d0*eta*c1y/den**(c1y+1.0d0)&
        + 4.0d0*eta*c1y/den**(c1y+1.0d0)&
        + 4.0d0*eta*eta2*c1y*(c1y+1.0d0)/den**(c1y+2.0d0)
    der2 = der2/(yend**2)
    detady(j) = 1.0d0/der1
    deta2dy2(j) = -der2/der1**3.0d0
END DO

DO j = 0, ny
    eta  = dfloat(j+ddrange(5))*dy + half*dy
    eta2 = (eta/yend)**2
    den = 2.0d0 - eta2
    ya(j) = eta/den**c1y
    der1 = 1.0d0/(den**c1y)+2.0d0*eta2*c1y/den**(c1y+1.0d0)
    der2 = 2.0d0*eta*c1y/den**(c1y+1.0d0)&
        + 4.0d0*eta*c1y/den**(c1y+1.0d0)&
        + 4.0d0*eta*eta2*c1y*(c1y+1.0d0)/den**(c1y+2.0d0)
    der2 = der2/(yend**2)
    detadya(j) = 1.0d0/der1
END DO
!
current_mesh => finest_mesh
dyf = current_mesh%dy
do j = -1, ny-1
    jj=j+1
    current_mesh%ey(j)=detady(jj)
    current_mesh%eya(j)=detadya(jj)
enddo

DO
    IF ( .NOT. ASSOCIATED(current_mesh%coarse_mesh) ) EXIT
    current_mesh => current_mesh%coarse_mesh

    ye = current_mesh%ye
    dyc = current_mesh%dy
    mny = anint(dyc/dyf)

    do j=-1,ye
        jj=j*mny+1
        if(jj > 0) then
            current_mesh%ey(j) = detady(jj)
            current_mesh%eya(j) = detady(jj+mny/2)
        else
            current_mesh%ey(j) = detady(0)
            current_mesh%eya(j) = detady(0)
        endif
         !!current_mesh%eya(j) = dyc/(y(jj+mny)-y(jj))
         !!current_mesh%ey(j) = (two*dyc)/(y(jj+mny)-y(jj-mny))
         !!current_mesh%eyy(j) = deta2dy2(jj)
    enddo
      !current_mesh%ey(0) =  detady(1)
      !current_mesh%eya(0) = detady(1+mny/2)
      !!current_mesh%eya(0) = dyc/(y(1+mny)-y(1))
      !!current_mesh%ey(0) =  dyc/(y(1+mny)-y(1))
      !!current_mesh%eyy(0) =  deta2dy2(1)

END DO

!-------------------------------------------------------------------
RETURN
END SUBROUTINE GRIDSETUP
!*******************************************************************

