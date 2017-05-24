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
!     Pressure Calculation
!
!     Note that the pressure only appears in the source
!       terms in the gas phase species and temperature
!       equations.
!
!     Input:  time - t
!             omg_p - frequency of pressure disturbance
!             epi_p - amplitude of pressure disturbance
!             period_p - period when no forcing is present
!
! ********************************************************************
REAL*8 FUNCTION PRESSURE(t)
  USE GLOBAL_DATA 
  USE MYMPI
  IMPLICIT NONE
!     Local Variables
  REAL*8  :: t
  
  PRESSURE = press
  
  RETURN
END FUNCTION PRESSURE
                         
REAL*8 FUNCTION DENSGAS(fvec,pscal)
  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE
  !Local Variables                                  
  REAL*8,intent(IN)  :: fvec(maxsize-1),pscal
  
  DENSGAS = dim_fact*pscal/fvec(1)

  RETURN
END FUNCTION DENSGAS
