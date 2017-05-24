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
!     Initial Conditions
!
!     This subroutine computes the initial conditions
!
! ********************************************************************
!
SUBROUTINE INITIAL_CONDITIONS1

USE GLOBAL_DATA
IMPLICIT NONE

!     Local Variables
INTEGER ::  i, k
      

!     Compute initial surface location
do k=drange(3)-1,drange(4)+1
    do i=drange(1)-1,drange(2)+1
        phi(i,k) = period_begin*period
        oldphi(i,k) = phi(i,k)
        dphidx(i,k) = 0.0d0
        dphidxa(i,k) = 0.0d0
        dphi2dx2(i,k) = 0.0d0
        dphidz(i,k) = 0.0d0
        dphidza(i,k) = 0.0d0
        dphi2dz2(i,k) = 0.0d0
    enddo
enddo
vec = 1.0d0
vvel = 0.5d0

return 
end SUBROUTINE INITIAL_CONDITIONS1
!*************************************************************

SUBROUTINE INITIAL_CONDITIONS2

USE GLOBAL_DATA
IMPLICIT NONE

!---------------------------------------------------------------
! Local Variables
INTEGER ::  i, k, j, eqn, col(7), ios, itemp
REAL*8  :: xi1, tinf, tsur, fact1, qsur, pressure
REAL*8 :: poo,coe,rhogas,t1,s1,s2,s3,t2
REAL*8 :: dpm(4),dpold(4)
CHARACTER*80 filnam

!  INITIALIZE (1) THE TEMPERATURE AND SPECIES CONCENTRATION
!             (2) THE BURN RATE
!             (3) THE MOMENTUM TERMS (VELOCITY)
!---------------------------------------------------------------

!!
dpold(1) = dble(ny)
dpold(2) = c1y
dpold(3) = yend
dpold(4) = press
itemp = 6
!!
! Interior grid
tinf = 2500.0d0
tsur = 810.0d0   !800
do j = drange(5), drange(6)
    do k = drange(3), drange(4)
        do i = drange(1), drange(2)
            xi1 = y(j)/0.01d0
            fact1 = dexp(-xi1)
            f(i,k,j,1) = tinf - (tinf - tsur)*fact1
            do eqn=2,neqmax-1
                f(i,k,j,eqn) = max( 0.4d0*xsurf(i,k,eqn)*fact1, 1.d-7)
            enddo
            f(i,k,j,neqmax) = tcold - (tcold - tsur)*fact1
        end do
    end do
end do
!!!

      col(1) = 1; col(2) = 2
      col(3)=1;col(4)=neqmax;
      col(5:6)=drange(5:6)
      col(7) = 0
      finest_mesh%w => f
      var_updated = 0
      call parallel_swap4(col,finest_mesh) 

      CALL  DIRICHLET_BC(itemp,finest_mesh)
      oldsoln = f
!
!     Compute initial burn rates
!
!!      CALL SMOOTH
      CALL SOLID_PROP
      CALL PYROLYSIS
!
!     Compute Initial Momemtntum
!
      poo = PRESSURE(tcyc)
      coe = dim_fact*poo
      q = zero
      if(ANY(type_bc(1,1,:) == 0) )then
         do j =0,nyv(1)+1
            if(type_bc(1,1,j) == 0) then
               do k = drange(3)-1, drange(4)+1
                  do i = -1 ,nxv(1)
                     fact1 = (f(i,k,j,1)+f(i+1,k,j,1))/&
                          &  (f(0,k,j,1)+f(1,k,j,1))
                     q(i,k,j,1) = vec_dirichlet(1,1,j)*fact1
                  enddo
               enddo
            endif
         enddo
      endif

      pold=zero; rate=zero
      do k = drange(3)-1, drange(4)+1
         do i = drange(1)-1, nxv(3)+1
            do j = drange(5), drange(6)
               rhogas = coe/(f(i,k,j,1))
               rhogas = two*coe/(f(i,k,j,1)+f(i,k,j+1,1))
               qsur = (rhos(i,k,0)/rhogas - one)*rb(i,k)
               q(i,k,j,3) = qsur
            end do
         end do
      end do


      col(1) = 1; col(2) = 2
      col(3)=1;col(4)=ndim;
      col(5:6)=drange(5:6)
      col(7) = 1
      finest_mesh%w => q
      var_updated = 1
      call parallel_swap4(col,finest_mesh) 
      do i = 1,ndim
         CALL  DIRICHLET_BC(i,finest_mesh)
      enddo


!----------------------------------------------------------------------
      RETURN 
    END SUBROUTINE INITIAL_CONDITIONS2
! ********************************************************************
