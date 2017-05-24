! ********************************************************************
! ********************************************************************
! Updated: TLJ; 1/12/2017
! Filename: ic.f90
! ********************************************************************
!
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
!
!     This subroutine computes the initial conditions
!
!     subroutines: INIT_FROM_SCRATCH
!                  INITIAL_CONDITIONS_PHI
!                  INITIAL_CONDITIONS_F
!                  INITIAL_CONDITIONS_Q
!
! ********************************************************************

SUBROUTINE INIT_FROM_SCRATCH

  USE GLOBAL_DATA
  IMPLICIT NONE

!---------------------------------------------------------------------
  INTEGER :: i,k,j
  REAL *8 :: burn_rate, coe,poo,qsur,pressure,rhogas
!---------------------------------------------------------------------

  ncyc = 0
  tinit = 0.0d0
  phi_begin = 0.0d0
  t_per = tinit
  phi_per = phi_begin

  if(ipack >= 0) CALL QUADRANTS
  CALL INITIAL_CONDITIONS_PHI
  if(ipack >= 0) THEN
     CALL  UPDATE_PSI(zero,zero)
  ELSE
     psi = -1.0d0     
  ENDIF
  CALL SOLID_PROP
  CALL PYROLYSIS   

  CALL INITIAL_CONDITIONS_F
  CALL INITIAL_CONDITIONS_Q

  CALL PYROLYSIS


!---------------------------------------------------------------------
  RETURN 
END SUBROUTINE INIT_FROM_SCRATCH
!*************************************************************
!
SUBROUTINE INITIAL_CONDITIONS_PHI

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

  call UPDATE_PSI(0,0)

  return 
end SUBROUTINE INITIAL_CONDITIONS_PHI
!*************************************************************

SUBROUTINE INITIAL_CONDITIONS_F

  USE GLOBAL_DATA
  IMPLICIT NONE

!---------------------------------------------------------------
! Local Variables
  INTEGER ::  i, k, j, eqn,l,jj,error
  REAL*8  :: xi1, tinf, tsur, fact1, qsur, pressure, poo,coe,rhogas,yy,tmp(10)

!  INITIALIZE (1) THE TEMPERATURE AND SPECIES CONCENTRATION
!             (2) THE BURN RATE
!             (3) THE MOMENTUM TERMS (VELOCITY)
!---------------------------------------------------------------

  ICDEBUG = .false.
  finest_mesh%w => f
  !call FourierExtVelUpdate(0,1)

!     Interior grid
  tinf = 3000.0d0
  tsur =  950.0d0
  !if (ipack==1 .and. NORADIATION.eqv..TRUE.) then
  !    tinf = 3000.0d0
  !    tinf = 2500.0d0
  !    tsur = 910.0d0
  !endif
!
  f(:,:,0,1) = tsur
  CALL SOLID_PROP
!
  do j = drange(5), drange(6)
     xi1 = y(j)/0.01d0
     fact1 = dexp(-xi1)
     do k = drange(3), drange(4)
        do i = drange(1), drange(2)
           f(i,k,j,1) = tinf - (tinf - tsur)*fact1
           do eqn=2,irxn
              f(i,k,j,eqn) = max(0.1d0*xsurf(i,k,eqn)*fact1,2.d-4)
           enddo
           do eqn = neqgas-1,neqgas
              f(i,k,j,eqn) = 1d-7
           end do

           f(i,k,j,neqmax) = tcold - (tcold - tsur)*fact1
        end do
     end do
  end do

  do eqn=1,neqmax
     CALL FILL_F_GHOST(eqn)
  enddo

  CALL SOLID_PROP
  CALL BC
!
! Compute initial burn rates
!
  CALL SOLID_PROP
  CALL PYROLYSIS
  CALL VELOCITY(0)
  call lambda(0)

! Recompute initial condition for T_solid
  if (NORADIATION.eqv..FALSE.) then
     do j = 0, ny
     do i = drange(1)-1,drange(2)+1
     do k = drange(3)-1,drange(4)+1
        f(i,k,j,neqmax) = tcold + ( f(i,k,0,1) - tcold) * &
             &exp(-vvel(i,k,0)*cp*y(j)/lambdas(i,k,j))
     enddo
     enddo
     enddo
  endif

  if (myid == 0) then
     write(*,'(80("&")/,a,/,20("-"))')'output from INITIAL_CONDITIONS_F:'
     write(6,*)'Initialization OK'
     write(*,'(20("-"),/a/,80("&"))')'END IC'
  endif


!----------------------------------------------------------------------
  RETURN 
END SUBROUTINE INITIAL_CONDITIONS_F
! ********************************************************************


! ********************************************************************
SUBROUTINE INITIAL_CONDITIONS_Q

  USE GLOBAL_DATA
  IMPLICIT NONE

!---------------------------------------------------------------
! Local Variables
  INTEGER ::  i, k, j, eqn
  REAL*8  :: qsur, pressure, poo,coe,rhogas
!---------------------------------------------------------------
!
! Compute Initial velocities
!
  poo = PRESSURE(tcyc(2))
  coe = one*dim_fact*poo

  q = zero
  pold = zero
  rate = zero
  do k = drange(3)-1, drange(4)+1
  do i = drange(1)-1, drange(2)+1
  do j = drange(5), drange(6)
     rhogas = coe/f(i,k,j,1)
     qsur = (rhos(i,k,0)/rhogas-1.0d0)*rb(i,k)
     q(i,k,j,3) = qsur
  end do
  end do
  end do

!----------------------------------------------------------------------
  RETURN 
END SUBROUTINE INITIAL_CONDITIONS_Q
! ********************************************************************
