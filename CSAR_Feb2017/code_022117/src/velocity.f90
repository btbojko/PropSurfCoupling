! ********************************************************************
! ********************************************************************
! Updated: TLJ; 12/29/2016
! Filename: velocity.f90
! ********************************************************************
! ********************************************************************
!
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
!
!     Subroutine computes the conserved velocities uvel (rho_gas*u),
!     vvel (rho_gas*vbar), wvel (rho_gas*w) from velocity vector q
!
!     subroutines: VELOCITY
!
!     Note: if Oseen, then vflag = 0
!
! ********************************************************************
SUBROUTINE VELOCITY(vflag)

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!---------------------------------------------------------------------
  INTEGER, INTENT(IN) :: vflag
! Local Variables
  INTEGER :: i, k, j, xs,zs,xe,ze,ys,ye,imax,kmax,jmax
  REAL*8  :: pressure,poo,coe,uu,vv,ww,ft, rhogas
  REAL*8  :: maxq(4),allmaxq(4),psimax,rbmax
  REAL*8  :: phi_x_avg,phi_z_avg
  logical outcome
!---------------------------------------------------------------------

  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  ys = drange(5)+1  !!vvel(0) is evaluated in bcq.f90
  ye = drange(6)

  poo = PRESSURE(tcyc(2))
  coe = two*dim_fact*poo

  mcomp = merge(2,1,vflag == 3)

  if (vflag.eq.0) then      !average at P points

     coe = dim_fact*poo
     do k=zs,ze
     do i=xs,xe
        do j=ys,ye
           MQterm(i,k,0) = (one +  MQchipr(j,1)*MQphi(i,k))
           if(irxn <=3) then
                 rhogas = coe/f(i,k,j,1)  ! <== this is probably rho = P/RT
           else
                 rhogas = rhog(i,k,j)
           end if
           uu = half*(q(i,k,j,1)+q(i-1,k,j,1))  ! these are cell center averages
           ww = half*(q(i,k,j,2)+q(i,k-1,j,2))
           vv = half*(q(i,k,j,3)+q(i,k,j-1,3))
           uvel(i,k,j) = rhogas*uu
           wvel(i,k,j) = rhogas*ww
           vvel(i,k,j) = rhogas*(MQvelocity(i,k,j,0)+vv-uu*MQchi(j,mcomp)*dphidx(i,k)&
                   &-ww*MQchi(j,mcomp)*dphidz(i,k))/MQterm(i,k,0)  ! see paper for v_bar
           vcnv(i,k,j) =  rhos(i,k,0)*rb(i,k)*sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)/rhogas 
           dtx(i,k,j) = dt/rhogas * time_coeff
           dtorhocp(i,k,j) = dt/rhogas/cp
        enddo
        ! burn surface kinematic equation
        ft0(i,k) = rb(i,k)*sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)   
        d_ft0(i,k) = d_rb(i,k)*sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)

        ! mass flux from surface
        vvel(i,k,0) = rhos(i,k,0)*ft0(i,k)
        d_vvel(i,k,0) = rhos(i,k,0)*d_ft0(i,k)
     enddo
     enddo

!--- OSEEN APPROX
     if(.NOT. doprojection) then
        !uvel(:,:,1:ny)=zero
        wvel(:,:,1:ny)=zero
        do k=zs,ze
        do i=xs,xe
        do j = drange(5),drange(6)  ! update over entire grid
           if(irxn <=3) then
              rhogas = coe/f(i,k,j,1)  ! <== this is probably rho = P/RT
           else
              rhogas = rhog(i,k,j)
           end if
           !>>> vic: shear, from juPack
           uvel(i, k, j) = rhogas*shear(j)
           !<<< vic
           MQterm(i,k,0) = (one +  MQchipr(j,1)*MQphi(i,k))
           vvel(i,k,j) = vvel(i,k,0)/MQterm(i,k,0)
           d_vvel(i,k,j) = d_vvel(i,k,0)/MQterm(i,k,0)
        enddo
        enddo
        enddo
     endif


  elseif (vflag.eq.1) then      !average at (u) points

     do j=ys,ye
     do k=zs,ze
     do i=xs,xe
        MQterm(i,k,1) = (one +  MQchipr(j,mcomp)*half*(MQphi(i,k)+MQphi(i+1,k)))

        rhogas = half*(rhog(i,k,j)+rhog(i+1,k,j)) 
        uu =  q(i,k,j,1)
        ww = quarter*(q(i,k,j,2) + q(i,k-1,j,2)+&
                   q(i+1,k,j,2)+q(i+1,k-1,j,2))
        vv = quarter*(q(i,k,j,3) + q(i,k,j-1,3)+&
                   q(i+1,k,j,3)+q(i+1,k,j-1,3))

        phi_z_avg =   half*(dphidz(i,k) + dphidz(i+1,k))

        vvel(i,k,j) = rhogas*(MQvelocity(i,k,j,1)+vv-uu*dphidxa(i,k)*MQchi(j,mcomp)&
                   &-ww*phi_z_avg*MQchi(j,mcomp))/MQterm(i,k,1)
        uvel(i,k,j) = rhogas*uu
        wvel(i,k,j) = rhogas*ww
!
!    dtx is dt/rho*time_coeff, since rho is averaged at different points
!    on the staggered grid it is reevaluated for each component
!
        dtx(i,k,j) = dt/rhogas*time_coeff

     enddo
     enddo
     enddo

  elseif (vflag.eq.2) then      !average at (w) points

     do j=ys,ye
     do k=zs,ze
     do i=xs,xe
        MQterm(i,k,2) = (one +  MQchipr(j,mcomp)*half*(MQphi(i,k)+MQphi(i,k+1)))

        rhogas = half*(rhog(i,k,j)+rhog(i,k+1,j)) 
        uu = quarter*(q(i,k,j,1) + q(i-1,k,j,1)+&
                   q(i,k+1,j,1)+q(i-1,k+1,j,1))
        ww =  q(i,k,j,2)
        vv = quarter*(q(i,k,j,3) + q(i,k,j-1,3)+&
                   q(i,k+1,j,3)+q(i,k+1,j-1,3))
        phi_x_avg =   half*(dphidx(i,k) + dphidx(i,k+1))

        vvel(i,k,j) = rhogas*(MQvelocity(i,k,j,2)+vv-uu*phi_x_avg*MQchi(j,mcomp)&
                   &-ww*dphidza(i,k)*MQchi(j,mcomp))/MQterm(i,k,2)
        uvel(i,k,j) = rhogas*uu
        wvel(i,k,j) = rhogas*ww
!
!    dtx is dt/rho*time_coeff, since rho is averaged at different points
!    on the staggered grid it is reevaluated for each component
!
        dtx(i,k,j) = dt/rhogas*time_coeff

     enddo
     enddo
     enddo

  elseif (vflag.eq.3) then      !average at (v) points

     do j=ys,ny-1
     do k=zs,ze
     do i=xs,xe 
        MQterm(i,k,3) = (one +  MQchipr(j,mcomp)*MQphi(i,k))

        rhogas = half*(rhog(i,k,j)+rhog(i,k,j+1)) 
        uu = quarter*(q(i,k,j,1) + q(i-1,k,j,1)+&
                   q(i,k,j+1,1)+q(i-1,k,j+1,1))
        ww = quarter*(q(i,k,j,2) + q(i,k-1,j,2)+&
                   q(i,k,j+1,2)+q(i,k-1,j+1,2))
        vv =  q(i,k,j,3)

        uvel(i,k,j) = rhogas*uu
        wvel(i,k,j) = rhogas*ww
        vvel(i,k,j) = rhogas*(MQvelocity(i,k,j,3)+vv-uu*dphidx(i,k)*MQchi(j,mcomp)&
                   &-ww*dphidz(i,k)*MQchi(j,mcomp))/MQterm(i,k,3)
!
!    dtx is dt/rho/time_coeff, since rho is averaged at different points
!    on the staggered grid it is reevaluated for each component
!
        dtx(i,k,j) = dt/rhogas*time_coeff

     enddo
     enddo
     enddo

  elseif (vflag.eq.4) then 
!
! PPE consrvetive average, drubar/dx+drvbar/dy+drwbar/dz = -dr/dt
!
     do k=zs,ze
     do i=xs,xe
        vvel(i,k,0) = rhos(i,k,0)*rb(i,k)&
                *sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)
     enddo
     enddo

     DO j=ys,ye

        do k = zs-1, ze+1
        do i = xs-1, xe+1
           MQterm(i,k,1) = (one +  MQchipr(j,1)*(MQphi(i+1,k)+MQphi(i,k))*half)
           MQterm(i,k,2) = (one +  MQchipr(j,1)*(MQphi(i,k+1)+MQphi(i,k))*half)
        enddo
        enddo
        do k=zs,ze
        do i=xs,xe

           quartuple =(/myid,i,k,j/)

           if(skipVarDensity) then
                 rhogas = rhoRef
           else
                 rhogas = coe/(f(i,k,j+1,1)+f(i,k,j,1))
                 if(j == ny) rhogas = coe/(two*f(i,k,j,1))
           endif
           uu = quarter*(q(i,k,j,1) + q(i-1,k,j,1)+&
                   q(i,k,j+1,1)+q(i-1,k,j+1,1))*MQchi(j,2)
           ww = quarter*(q(i,k,j,2) + q(i,k-1,j,2)+&
                   q(i,k,j+1,2)+q(i,k-1,j+1,2))*MQchi(j,2)

           rate(i,k,j,1) = quarter*(q(i,k,j,1) + q(i-1,k,j,1)+&
                   q(i,k,j+1,1)+q(i-1,k,j+1,1))
           rate(i,k,j,2) = quarter*(q(i,k,j,2) + q(i,k-1,j,2)+&
                   q(i,k,j+1,2)+q(i,k-1,j+1,2))

           vv =  q(i,k,j,3)
           vvel(i,k,j) = rhogas*(MQvelocity(i,k,j,3)+vv-uu*dphidx(i,k)-ww*dphidz(i,k))

        enddo
        enddo

        do k=zs-1,ze
        do i=xs-1,xe
           quartuple =(/myid,i,k+1,j/)
           if(skipVarDensity) then
                 rhogas = rhoRef
           else
                 rhogas = coe/(f(i,k,j,1)+f(i+1,k,j,1)) 
           endif
           uvel(i,k,j) = rhogas*q(i,k,j,1)*MQterm(i,k,1)

           if(skipVarDensity) then
                 rhogas = rhoRef
           else
                 rhogas = coe/(f(i,k,j,1)+f(i,k+1,j,1)) 
           endif
           wvel(i,k,j) = rhogas*q(i,k,j,2)*MQterm(i,k,2)
        enddo
        enddo

     ENDDO

  elseif (vflag.eq.-4) then 
!
! PPE consrvetive average, drubar/dx+drvbar/dy+drwbar/dz = -dr/dt
!
     DO j=ys,ye

        do k=zs,ze
        do i=xs,xe

           if(skipVarDensity) then
                 rhogas = rhoRef
           else
                 rhogas = coe/(f(i,k,j+1,1)+f(i,k,j,1))
                 if(j == ny) rhogas = coe/(two*f(i,k,j,1))
           endif
           uu = quarter*(q(i,k,j,1) + q(i-1,k,j,1)+&
                   q(i,k,j+1,1)+q(i-1,k,j+1,1))
           ww = quarter*(q(i,k,j,2) + q(i,k-1,j,2)+&
                   q(i,k,j+1,2)+q(i,k-1,j+1,2))
           vvel(i,k,j) = rhogas*(uu*dphidx(i,k)+ww*dphidz(i,k))

        enddo
        enddo

     ENDDO

  elseif (vflag.eq.5) then
!TEST

     DO j=ys,ye

        do k = zs-1, ze+1
        do i = xs-1, xe+1
           MQterm(i,k,1) = (one +  MQchipr(j,1)*(MQphi(i+1,k)+MQphi(i,k))*half)
           MQterm(i,k,2) = (one +  MQchipr(j,1)*(MQphi(i,k+1)+MQphi(i,k))*half)
        enddo
        enddo
        do k=zs,ze
        do i=xs,xe

           if(skipVarDensity) then
                 rhogas = rhoRef
           else
                 rhogas = coe/(f(i,k,j+1,1)+f(i,k,j,1))
                 if(j == ny) rhogas = coe/(two*f(i,k,j,1))
           endif
           uu = rate(i,k,j,1)*MQchi(j,2)

           ww = rate(i,k,j,2)*MQchi(j,2)
           vv = q(i,k,j,3)
           vvel(i,k,j) = rhogas*(MQvelocity(i,k,j,3)+vv-uu*dphidx(i,k)-ww*dphidz(i,k))

        enddo
        enddo

        do k=zs-1,ze
        do i=xs-1,xe
           quartuple =(/myid,i,k+1,j/)
           if(skipVarDensity) then
                 rhogas = rhoRef
           else
                 rhogas = coe/(f(i,k,j,1)+f(i+1,k,j,1)) 
           endif
           uvel(i,k,j) = rhogas*q(i,k,j,1)*MQterm(i,k,1)

           if(skipVarDensity) then
                 rhogas = rhoRef
           else
                 rhogas = coe/(f(i,k,j,1)+f(i,k+1,j,1)) 
           endif
           wvel(i,k,j) = rhogas*q(i,k,j,2)*MQterm(i,k,2)

        enddo
        enddo

     ENDDO

  elseif (vflag.eq.6) then 
!
! TEST ONLY
!
     DO j=ys,ye

        do k=zs,ze
        do i=xs,xe

           rhogas = coe/(f(i,k,j+1,1)+f(i,k,j,1))
           if(j == ny) rhogas = coe/(two*f(i,k,j,1))
           uu = uvel(i,k,j)
           ww = wvel(i,k,j)
           vv = q(i,k,j,3)
           vvel(i,k,j) = rhogas*(MQvelocity(i,k,j,3)+vv-uu*dphidx(i,k)-ww*dphidz(i,k))

        enddo
        enddo

        do k=zs-1,ze
        do i=xs-1,xe

           rhogas = coe/(f(i,k,j,1)+f(i+1,k,j,1)) 
           uvel(i,k,j) = rhogas*q(i,k,j,1)
           rhogas = coe/(f(i,k,j,1)+f(i,k+1,j,1)) 
           wvel(i,k,j) = rhogas*q(i,k,j,2)

        enddo
        enddo
!
     ENDDO

  elseif (vflag.eq.7) then
!
!   VISUALIZATION ONLY, vflag ==7 is used for visualization purposes to put
!   the vectors back in the cell centered location so to visualize them using
!   either MATLAB OR ROCKETEER
!

     maxq(1) = 1000.
     do k=zs,ze
     do i=xs,xe
     do j=ys,ye

        uu = half*(q(i,k,j,1)+q(i-1,k,j,1))
        ww = half*(q(i,k,j,2)+q(i,k-1,j,2))
        vv = half*(q(i,k,j,3)+q(i,k,j-1,3))
        uvel(i,k,j) = uu
        wvel(i,k,j) = ww
        vvel(i,k,j) = vv

     enddo
     enddo
     enddo

  endif  !vflag.eq.0

  RETURN
END SUBROUTINE VELOCITY
!**************************************************************************
