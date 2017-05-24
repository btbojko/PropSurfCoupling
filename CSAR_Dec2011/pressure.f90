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
!--------------------------------------------
  USE GLOBAL_DATA 
  USE MYMPI
  IMPLICIT NONE
  
!     Local Variables
  REAL*8  ::  t, per_test, delta_p
!--------------------------------------------

  phi0 = phi(0,0)
  call MPI_BCAST(phi0,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)

!  pressure = press
!  pressure = press*(1.0d0+epi_p*sin(2.0d0*pi*omg_p*t/period_p))
  pressure = press*( 1.0d0 + epi_p * sin(2.0d0*pi*omg_p*t &
       + phase_p*half*pi) )

  if(tcyc < time_steady_pressure) pressure = press

  dpress_dt =  zero

!  if (ipack.eq.1) then
!     per_test = dint(phi0/period)
!     if (per_test .ge. nperiod-2) then
!        if (iper.eq.0) then
!           period_p = t
!           iper = 1
!           write(6,*)'period_p,iper=', period_p,iper
!        endif
!     endif
!     if (per_test .ge. nperiod-1) then
!        if (iper.eq.1) then
!           period_p = 0.5*(period_p + 2.0*(t-period_p)+t)
!           iper = 2
!           write(6,*)'period_p,iper=', period_p,iper
!        endif
!        pressure = 20.0 + 20.0*exp(-(t-period_p)**2/0.00001334)
!        pressure = 0.5*(60.0 + 20.0*tanh((t-period_p)/0.001))
!        pressure = 0.5*(60.0 - 20.0*tanh((t-period_p)/0.001))
!        pressure = 20.0 + 20.0*exp(-(t-period_p)**2/delta_p)
!        dpress_dt =  -20.0*2.0d0*(t - period_p)/delta_p*&
!             exp(-(t-period_p)**2/delta_p)
!     endif
!  endif

  RETURN
END FUNCTION PRESSURE
