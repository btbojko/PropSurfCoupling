! ********************************************************************
! ********************************************************************
! Updated: TLJ; 12/19/2016
! Filename: conductivity.f90
! ********************************************************************
! ********************************************************************
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
!
!     This set of subroutines computes the gas-phase thermal
!     conductivity lambda_gas, as well as the volume fraction
!     of aluminum (if present).
!
!     subroutines: LAMBDA
!                  AL_VF
!
! ********************************************************************

SUBROUTINE LAMBDA(lavg)

  USE GLOBAL_DATA

  IMPLICIT NONE
!----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: lavg
! Local Variables
  INTEGER :: eqn,i,j,k,ys,ye
  REAL*8  :: hx, hy, hz
  REAL*8  :: coe, poo, tavg(3),pressure
!----------------------------------------------------------------------

!
! now update gas-phase lambda
!
  poo = PRESSURE(tcyc(2))
  coe = two*dim_fact*poo

  hx = 1.0d0/(2.0d0*dx)
  hy = 1.0d0/(2.0d0*dy)
  hz = 1.0d0/(2.0d0*dz)

  if (lavg.eq.0) then

     do j=drange(5),drange(6)
     do k=drange(3)-1,drange(4)+1
     do i=drange(1)-1,drange(2)+1

        lambdag(i,k,j) = a_lamg*f(i,k,j,1)**e_lamg + b_lamg
        rhog(i,k,j)    =  coe/f(i,k,j,1)

     enddo
     enddo
     enddo

     if (NORADIATION.eqv..FALSE.) call Al_VF

     ys = drange(5)
     ye = drange(6)

     if(dSerialRange(5).eq.0) ys = ys + 1
     if(dSerialRange(6).eq.ny) ye = ye - 1

!
! Now compute derivatives for interior grid
!
     do j = ys, ye
     do k = drange(3), drange(4)
     do i = drange(1), drange(2)
        dlgdx(i,k,j) = (lambdag(i+1,k,j)-lambdag(i-1,k,j))*hx
        dlgdz(i,k,j) = (lambdag(i,k+1,j)-lambdag(i,k-1,j))*hz
        dlgdy(i,k,j) = (lambdag(i,k,j+1)-lambdag(i,k,j-1))*hy
     end do
     end do
     end do

  ELSEif (lavg.eq.1) then

     do j=drange(5),drange(6)
     do k=drange(3)-1,drange(4)+1
     do i=drange(1)-1,drange(2)+1

        lambdag(i,k,j) = a_lamg*tmdpnt(i,k,j)**e_lamg + b_lamg

     enddo
     enddo
     enddo

     do i=drange(1)-1,drange(2)
     do k=drange(3)-1,drange(4)
     do j=drange(5)-1,drange(6)

        tavg(1) = quarter*(tmdpnt(i,k,j) + tmdpnt(i+1,k,j)+&   ! i+1/2,  k+1/2  ,j
                   tmdpnt(i,k+1,j)+tmdpnt(i+1,k+1,j))
        tavg(2) = quarter*(tmdpnt(i,k,j) + tmdpnt(i+1,k,j)+&   ! i+1/2,  k      ,j+1/2
                   tmdpnt(i,k,j+1)+tmdpnt(i+1,k,j+1))
        tavg(3) = quarter*(tmdpnt(i,k,j) + tmdpnt(i,k+1,j)+&   ! i,     k+1/2   ,j+1/2
                   tmdpnt(i,k,j+1)+tmdpnt(i,k+1,j+1))

        lbdavg(i,k,j,:) =  a_lamg*tavg(:)**e_lamg + b_lamg

     enddo
     enddo
     enddo


  endif

!---------------------------------------------------------------------
  RETURN
END SUBROUTINE LAMBDA
!*********************************************************************
!*********************************************************************
subroutine Al_VF

  use global_data
  USE MYMPI
  implicit none

  integer :: i,k,j,m
  real*8,Dimension(maxsize) :: Wfc,rhc,Vfc,lmc
  real*8:: sumvol
!--------------------------------------------------

  do j = drange(5),drange(6)
  do k = drange(3)-1,drange(4)+1
  do i = drange(1)-1,drange(2)+1
     do m = 2,3
              Wfc(m) = max(f(i,k,j,3+m),1d-20)
     end do
     if(sum(Wfc(2:3)) > 1d0) then
              Wfc(2:3) = Wfc(2:3) + (1d0-sum(Wfc(2:3)))*0.55d0
     end if
     Wfc(1) = 1d0-sum(Wfc(2:3))

     rhc(1:3) = (/rhog(i,k,j),rho_al,rho_ala/)
     Vfc(1:3) = Wfc(1:3)/rhc(1:3)
     sumvol = sum(Vfc(1:3))
     Vfrac(i,k,j,1:2)  = Vfc(2:3)/sumvol
     rhog(i,k,j) = sum(Wfc(1:3))/sumvol
     lmc(1:3) = (/lambdag(i,k,j),lambda_al,lambda_ala/)
     !circular shifts
     Vfc(1:3) = cshift(Vfc(1:3),-1)
     Lmc(1:3) = cshift(Lmc(1:3),-1)
     Rhc(1:3) = cshift(Rhc(1:3),-1)
     tmdpnt(i,k,j) = lambdag(i,k,j)
     !sumvol is a junk output
     CALL BLEND(Vfc(1:3),lmc(1:3),rhc(1:3),lambdag(i,k,j),sumvol,3)  
     DiamP(i,k,j,1) = initial_diameter*abs(f(i,k,j,neqgas-1)*&
          beta_al/(f(i,k,j,neqgas-1)*beta_al &
          +f(i,k,j,neqgas)))**(1d0/3d0) * 1d-4
     DiamP(i,k,j,2) = initial_diameter*abs(f(i,k,j,neqgas)/&
          (f(i,k,j,neqgas-1)*beta_al  +f(i,k,j,neqgas)) &
          * beta_al/beta_rhoal)**(1d0/3d0) * 1d-4
  end do
  end do
  end do

  do j = lbound(VOD,3),ubound(VOD,3)  
  do k = drange(3)-1,drange(4)+1
  do i = drange(1)-1,drange(2)+1
     do m = 1,2
        if(sum(f(i,k,j,neqgas-1:neqgas)) > 1d-12) then
           if(abs(DiamP(i,k,j,m)) > 1d-14) then
              Vod(i,k,j,m) = Vfrac(i,k,j,m)/DiamP(i,k,j,m)
           else
              Vod(i,k,j,m) = 0d0
           end if
        else
           Vod(i,k,j,m) = 0d0
        end if
     end do
  end do
  end do
  end do

  do k = drange(3)-1,drange(4)+1
  do i = drange(1)-1,drange(2)+1
     do m = 1,2
        Vfrac(i,k,lbound(Vfrac,3):drange(5),m) = Vfrac(i,k,drange(5),m)
        VOD(i,k,lbound(VOD,3):drange(5),m) = VOD(i,k,drange(5),m)
     end do
  end do
  end do
              
           
  call mpi_barrier(comm3d,i)  !cancel
! 
  return
end subroutine Al_VF
!*********************************************************************
