! ********************************************************************
! ********************************************************************
! Updated: TLJ; 2/1/2017
! Filename: flux.f90
! ********************************************************************
!
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
!
!     This subroutine updates the phi function using WENO
!       for the fluxes
!
!     subroutines: FLUX
!
! ********************************************************************

SUBROUTINE FLUX(t,tout)

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!---------------------------------------------------------------------
! Dummy Variables
  REAL*8, INTENT(IN) :: t,tout
! Local Variables
  INTEGER ::  i, k, nspeed, all_nspeed
  REAL*8  :: dtt,eps,mymax,mymin
  REAL*8  :: fx,fz,is1,is2,ap1,ap2,phix1,phix2,phiz1,phiz2
  REAL*8  :: all_speed
!---------------------------------------------------------------------

! compute time step
  dtt = dt

  eps = 1.0d-8
  do k=drange(3),drange(4)
  do i=drange(1),drange(2)
     !
     ! X-DERIVATIVE CALCULATION
     !
     fx = 0.d0
     if(phi(i-1,k).le.phi(i+1,k).and.&
             phi(i-1,k).lt.phi(i,k)) then
        is1 = (phi(i-2,k)-2.*phi(i-1,k)+phi(i,k)  )**2
        is2 = (phi(i-1,k)-2.*phi(i,k)  +phi(i+1,k))**2
        ap1 = 0.5/(eps+is1)
        ap2 = 1.0/(eps+is2)
        phix1 = (phi(i-2,k)-4.*phi(i-1,k)+3.*phi(i,k))/(2.*dx)
        phix2 = (phi(i+1,k)-phi(i-1,k))/(2.*dx)
        fx = (ap1*phix1+ap2*phix2)/(ap1+ap2)
     endif

     if(phi(i+1,k).le.phi(i-1,k).and.&
             phi(i+1,k).lt.phi(i,k)) then
        is1 = (phi(i-1,k)-2.*phi(i,k)  +phi(i+1,k))**2
        is2 = (phi(i,k)  -2.*phi(i+1,k)+phi(i+2,k))**2
        ap1 = 1.0/(eps+is1)
        ap2 = 0.5/(eps+is2)
        phix1 = (phi(i+1,k)-phi(i-1,k))/(2.*dx)
        phix2 = (-3.*phi(i,k)+4.*phi(i+1,k)-phi(i+2,k))/(2.*dx)
        fx = (ap1*phix1+ap2*phix2)/(ap1+ap2)
     endif
     !
     ! Z-DERIVATIVE CALCULATION
     !
     fz = 0.d0
     if(.not. is2d) then
        if(phi(i,k-1).le.phi(i,k+1).and.&
                phi(i,k-1).lt.phi(i,k) ) then
           is1 = (phi(i,k-2)-2.*phi(i,k-1)+phi(i,k)  )**2
           is2 = (phi(i,k-1)-2.*phi(i,k)  +phi(i,k+1))**2
           ap1 = 0.5/(eps+is1)
           ap2 = 1.0/(eps+is2)
           phiz1 = (phi(i,k-2)-4.*phi(i,k-1)+3.*phi(i,k))/(2.*dz)
           phiz2 = (phi(i,k+1)-phi(i,k-1))/(2.*dz)
           fz = (ap1*phiz1+ap2*phiz2)/(ap1+ap2)
        endif

        if(phi(i,k+1).le.phi(i,k-1).and.&
                phi(i,k+1).lt.phi(i,k) ) then
           is1 = (phi(i,k-1)-2.*phi(i,k)  +phi(i,k+1))**2
           is2 = (phi(i,k)  -2.*phi(i,k+1)+phi(i,k+2))**2
           ap1 = 1.0/(eps+is1)
           ap2 = 0.5/(eps+is2)
           phiz1 = (phi(i,k+1)-phi(i,k-1))/(2.*dz)
           phiz2 = (-3.*phi(i,k)+4.*phi(i,k+1)-phi(i,k+2))/(2.*dz)
           fz = (ap1*phiz1+ap2*phiz2)/(ap1+ap2)
        endif
     end if

     phit(i,k) = rb(i,k)*sqrt(one+fx*fx+fz*fz)
  enddo
  enddo
  !
  ! Update phi
  !
  do k=drange(3),drange(4)
  do i=drange(1),drange(2)
     phi(i,k) = phi(i,k) - dtt*phit(i,k)
  enddo
  enddo
  !
  ! EVALUATE THE AVERAGE SPEED
  !
  speed = 0
  nspeed = 0
  do i = drange(1),drange(2)
  do k = drange(3),drange(4)
     speed = speed + phit(i,k)
     nspeed = nspeed + 1
  enddo
  enddo

  call MPI_ALLREDUCE(speed,all_speed,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
  call MPI_ALLREDUCE(nspeed,all_nspeed,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
  speed = all_speed/all_nspeed

  speed_timeint = speed_timeint + speed * dtt
  dt_timeint = dt_timeint + dtt

!---------------------------------------------------------------------
  RETURN 
END SUBROUTINE FLUX
!*********************************************************************
