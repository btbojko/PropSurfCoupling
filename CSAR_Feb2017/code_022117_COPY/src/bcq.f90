! ********************************************************************
! ********************************************************************
! Updated: TLJ; 12/29/2016
! Filename: bcq.f90
! ********************************************************************
! ********************************************************************
!
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
!
!     Subroutine updates boundary conditions for velocity
!
!     subroutines: BCQ
!                  FourierExtVelUpdate; ignored at the moment
!
! ********************************************************************

SUBROUTINE BCQ(flag)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE
  include "parallel.h"

  !----------------------------------------------------------------
  ! Local variables
  INTEGER,INTENT(IN) :: flag
  INTEGER :: i, k, j, eqn,xs,zs,xe,ze,col(7)
  REAL *8 :: poo,coe,usurf,wsurf,vsurf
  REAL *8 :: pressure, rratio, my_flux, dya, rhogas
  REAL *8 :: cjswall,cj1wall,cj2wall
  !----------------------------------------------------------------

  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  poo = PRESSURE(tcyc(2))
  coe = two*dim_fact*poo
  col(1)=1; col(2)=2; col(3)=1; col(4) = ndim;
  col(5)=drange(5); col(6)=col(5)      !j=0 ONLY!
  col(7) = 0

  !
  !  IMPOSE velocity boundary conditions
  !

  ! Note: for Oseen, q is never used, and uvel(ny),wvel(ny)
  !  are also not used since the equations for T and species
  !  use zero gradients at j=ny, and uses jump conditions
  !  at j=0 which are again independent of velocities.

  do i=xs-1,xe+1
  do k=zs-1,ze+1

     !SOUTH; updates vvel only

     rhogas=coe/(f(i,k,0,1)+f(i,k,1,1)) !rhogas on the surf
     rratio = rhos(i,k,0)/rhogas-1.0d0
     vsurf=rb(i,k)**2*rratio/vec(i,k,0)    
     if(vsurf<0) write(*,*)i,k,'ERROR vsur <0 at wall',vsurf,&
         myid,'FVECT',f(i,k,0:1,1:5),rratio,rb(i,k),rhogas
     newqbnd(i,k,3) = vsurf
     vvel(i,k,0) = rhos(i,k,0)*vec(i,k,0)    !rhogas*vvel

     !NORTH

     vvel(i,k,ny) = vvel(i,k,ny-1)
     q(i,k,ny:ny+1,1:2) = 0.0d0
     uvel(i,k,ny) = 0.0d0
     wvel(i,k,ny) = 0.0d0

!........................................................................
! These formulae preserve the derivative at the boundary j=ny
!........................................................................
!     q(i,k,ny+1,1) = q(i,k,ny-2,1)-3.0*q(i,k,ny-1,1)+3.0*q(i,k,ny,1)
!     q(i,k,ny+1,2) = q(i,k,ny-2,2)-3.0*q(i,k,ny-1,2)+3.0*q(i,k,ny,2)
!
  enddo
  enddo

  ! South, sets [qxn]=0; not used if Oseen
  do i=xs-1,xe
  do k=zs-1,ze
     usurf = -half*(q(i+1,k,0,3)+q(i,k,0,3))*dphidxa(i,k)  !-vsurf*fx
     wsurf = -half*(q(i,k,0,3)+q(i,k+1,0,3))*dphidza(i,k)  !-vsurf*fz
     q(i,k,0,1:2) = (/usurf, wsurf/)
  enddo
  enddo

  !
  !  This extra communication can be reduced
  !
  finest_mesh%w => q
  call parallel_swap4(col,finest_mesh)

  cjswall=8.0d0/3.0d0
  cj1wall=two
  cj2wall=third
  do i=xs-1,xe+1
  do k=zs-1,ze+1
     uvel(i,k,0) = q(i,k,0,1)
     wvel(i,k,0) = q(i,k,0,2)
     srfqbnd(i,k,1) = uvel(i,k,0)
     srfqbnd(i,k,2) = wvel(i,k,0)
!
!   2nd order extrapolation at the wall, file BCUW.maple
!
     newqbnd(i,k,1) = cjswall*srfqbnd(i,k,1) - cj1wall*q(i,k,1,1) + cj2wall*q(i,k,2,1)
     newqbnd(i,k,2) = cjswall*srfqbnd(i,k,2) - cj1wall*q(i,k,1,2) + cj2wall*q(i,k,2,2)
  enddo
  enddo
  
  if(flag == 1)then   !reset values
     do i=xs-1,xe+1
        do k=zs-1,ze+1
           q(i,k,0,1:3) = oldsoln(i,k,0,maxsize+1:maxsize+3)
        enddo
     enddo
  elseif(flag == 2)then  !plug in new values
     do i=xs-1,xe+1
        do k=zs-1,ze+1
           q(i,k,0,1:3) =  newqbnd(i,k,1:3)
        enddo
     enddo
  endif

  !-------------------------------------------------------------------------------
  RETURN
END SUBROUTINE BCQ
!*******************************************************************
!*******************************************************************
SUBROUTINE FourierExtVelUpdate(flag,eqn)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE
  include "parallel.h"

  !----------------------------------------------------------------
  ! Local variables
  INTEGER :: j,iMed,flag,id,eqn
  REAL *8 :: YRFLO,distloc
  REAL *8 :: rhoC,stdImposed
  REAL *8 :: meanProcess,stdProcess,UnifRnd(2),TwoPi
  !----------------------------------------------------------------

  if(skipCrossflow) then
     YRFLO = 1d99
     jRFLO(1:3) = ny+1
     !no crossflow
     q(:,:,ny:ny+1,1:2) = zero
     q(:,:,ny,3) = q(:,:,ny-1,3)
     q(:,:,ny+1,3) = q(:,:,ny-1,3)
     if(allocated(Uimposed)) Uimposed = zero
     return
  endif
  TRFLO = 2.8d3
  !
  !
  YRFLO = 0.064942d0
  jRFLO = ny+1
  distloc = 1d99
  DO j = 0, ny+1
     if(abs(y(j) - YRFLO) < distloc) then
        distloc = abs(y(j) - YRFLO )
        jRFLO(1:3) = j
     endif
  ENDDO

  if(flag  == 0) then
     Uimposed(1:2) = matrix%m(1)*cos(tcyc(1:2)*2*pi*matrix%f(1) + matrix%a(1))
     do id = 2,matrix%nFourier
        Uimposed(1:2) = Uimposed(1:2) + matrix%m(id)*cos(tcyc(1:2)*2*pi*matrix%f(id) + matrix%a(id));
     enddo

     Uimposed = (Uimposed + Mimposed)*m2cm

  elseif(flag == 1) then

     !this should be done on the last iteration within a timestep only
     q(:,:,jRFLO(1):ny+1,1) = Uimposed(1)
     q(:,:,jRFLO(2):ny+1,2) = zero
     do j = jRFLO(3),ny+1
        q(:,:,j,3) = q(:,:,jRFLO(3)-1,3)
     enddo

  elseif(flag == 2) then

     if(eqn == 1)then
        q(:,:,jRFLO(1):ny+1,1) = Uimposed(1)  !OLD value
     elseif(eqn == 2)then
        q(:,:,jRFLO(2):ny+1,2) = zero
     elseif(eqn == 3)then
        do j = jRFLO(3),ny+1
           q(:,:,j,3) = q(:,:,jRFLO(3)-1,3)
        enddo
     endif

     dqdt(:,:,jRFLO(eqn):ny+1,eqn)     = 0d0
     dconv(:,:,jRFLO(eqn):ny+1,1:3)    = 0d0
     Bm(1:2,1:4,:,:,jRFLO(eqn):ny+1)     = 0d0

  elseif(flag == 3) then

     divt(:,:,jRFLO(1)+1:ny+1) = 0d0

  elseif(flag == 4) then

     p(:,:,jRFLO(1)-1:ny) = 0d0
     wrk(:,:,jRFLO(1)-1:ny) = 0d0

  elseif(flag == 5) then

     !this should be done on the last iteration within a timestep only
     q(:,:,jRFLO(1):ny+1,1) = Uimposed(2)
     q(:,:,jRFLO(2):ny+1,2) = zero
     do j = jRFLO(3),ny+1
        q(:,:,j,3) = q(:,:,jRFLO(3)-1,3)
     enddo
  endif

  !-------------------------------------------------------------------------------
  RETURN
END SUBROUTINE FourierExtVelUpdate
!*******************************************************************

