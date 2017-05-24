!**********************************************************************
SUBROUTINE VELOCITY(vflag)
  
USE GLOBAL_DATA
USE MYMPI
IMPLICIT NONE
  
!---------------------------------------------------------------------
INTEGER, INTENT(IN) :: vflag
! Local Variables
INTEGER :: i, ig, k, j, xs,zs,xe,ze,ys,ye,imax,kmax,jmax
REAL*8  :: pressure,poo,coe_rho,uu,vv,ww
REAL*8  :: ft,phi_x_avg,phi_z_avg,rhogas
!---------------------------------------------------------------------
  
xs = drange(1)
xe = drange(2)
zs = drange(3)
ze = drange(4)
ys = drange(5)+1  !!vvel(0) is evaluated in the boundary conditions
ye = nyv(vflag)+1

poo = PRESSURE(tcyc)
coe_rho = dim_fact*poo

if (vflag.eq.0) then      !average at P points

    do k=zs,ze
        do i=xs,nxv(0)
            do j=ys,nyv(0)
                rhogas = coe_rho/f(i,k,j,1)
                uu = half*(q(i,k,j,1)+q(i-1,k,j,1))
                ww = half*(q(i,k,j,2)+q(i,k-1,j,2))
                vv = half*(q(i,k,j,3)+q(i,k,j-1,3))
                uvel(i,k,j) = rhogas*uu
                wvel(i,k,j) = rhogas*ww
                vvel(i,k,j) = rhogas*(rb(i,k)*vec(i,k)+vv-uu*dphidx(i,k)-ww*dphidz(i,k))
                vcnv(i,k,j) =  rhos(i,k,0)*rb(i,k)&
                    *sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)/rhogas 
                dtx(i,k,j) = dt/rhogas*time_coeff
            enddo
            vvel(i,k,0) = rhos(i,k,0)*rb(i,k)&
                *sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)
            ft0(i,k) = rb(i,k)*sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)
            d_ft0(i,k) = d_rb(i,k)*sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)
        enddo
    enddo
 
!--- OSEEN APPROX
    if(.NOT. doprojection) then
        uvel(:,:,1:ny)=zero;wvel(:,:,1:ny)=zero
        do k=zs,ze
            do i=xs,xe
                vvel(i,k,:) = vvel(i,k,0)
            enddo
        enddo
    endif

elseif (vflag.eq.1) then      !average at (u) points

    do j=ys,nyv(vflag)+1
        do k=zs,ze
            do i=xs-1,nxv(vflag)+1
                rhogas = half*(rhog(i,k,j)+rhog(i+1,k,j)) 
                uu =  q(i,k,j,1)
                ww = quarter*(q(i,k,j,2) + q(i,k-1,j,2)+&
                    q(i+1,k,j,2)+q(i+1,k-1,j,2))
                vv = quarter*(q(i,k,j,3) + q(i,k,j-1,3)+&
                    q(i+1,k,j,3)+q(i+1,k,j-1,3))
                phi_z_avg =   half*(dphidz(i,k) + dphidz(i+1,k))
                ft = half*(rb(i,k)*vec(i,k)+rb(i+1,k)*vec(i+1,k))
                vvel(i,k,j) = rhogas*(ft+vv-uu*dphidxa(i,k)-ww*phi_z_avg)
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

    do j=ys,nyv(vflag)+1
        do k=zs,ze
            do i=xs-1,nxv(vflag)+1

                rhogas = half*(rhog(i,k,j)+rhog(i,k+1,j)) 
                uu = quarter*(q(i,k,j,1) + q(i-1,k,j,1)+&
                    q(i,k+1,j,1)+q(i-1,k+1,j,1))
                ww =  q(i,k,j,2)
                vv = quarter*(q(i,k,j,3) + q(i,k,j-1,3)+&
                    q(i,k+1,j,3)+q(i,k+1,j-1,3))
                ft = half*(rb(i,k)*vec(i,k)+rb(i,k+1)*vec(i,k+1))
                vvel(i,k,j) = rhogas*(ft+vv-uu*dphidx(i,k)-ww*dphidz(i,k))
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

    do j = ys, nyv(vflag)+1
        do k = zs-1, ze+1
            do i = xs-1,nxv(vflag)+1
                rhogas = half*(rhog(i,k,j)+rhog(i,k,j+1)) 
                uu = quarter*(q(i,k,j,1) + q(i-1,k,j,1)+&
                    q(i,k,j+1,1)+q(i-1,k,j+1,1))
                ww = quarter*(q(i,k,j,2) + q(i,k-1,j,2)+&
                    q(i,k,j+1,2)+q(i,k-1,j+1,2))
                vv =  q(i,k,j,3)
                uvel(i,k,j) = rhogas*uu
                wvel(i,k,j) = rhogas*ww
                vvel(i,k,j) = rhogas*(rb(i,k)*vec(i,k)+vv-uu*dphidx(i,k)-ww*dphidz(i,k))
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
    DO j=ys,ny
        do k = 0, finest_mesh%ze
            do i = 0, finest_mesh%xe+1
                rhogas = two*coe_rho/(f(i,k,j+1,1)+f(i,k,j,1))
                if(j == ny) rhogas = coe_rho/f(i,k,j,1)
                uu = quarter*(q(i,k,j,1) + q(i-1,k,j,1)+&
                    q(i,k,j+1,1)+q(i-1,k,j+1,1))
                ww = quarter*(q(i,k,j,2) + q(i,k-1,j,2)+&
                    q(i,k,j+1,2)+q(i,k-1,j+1,2))
                vv =  q(i,k,j,3)
                    vvel(i,k,j) = rhogas*(rb(i,k)*vec(i,k)+vv-uu*dphidx(i,k)-ww*dphidz(i,k))
            enddo
        enddo
        do k= -1,finest_mesh%ze+1
            do i= -1,finest_mesh%xe+1
                rhogas = two*coe_rho/(f(i,k,j,1)+f(i+1,k,j,1))
                uvel(i,k,j) = rhogas*q(i,k,j,1)
                rhogas = two*coe_rho/(f(i,k,j,1)+f(i,k+1,j,1)) 
                wvel(i,k,j) = rhogas*q(i,k,j,2)
            enddo
        enddo

        if(ALL(type_bc(3,3,:) == 0) .AND. neqmax == neqgas +1) then  !i'm saying there is solid
            do k=zs,ze
                do i=0,nxv(3)
                    vvel(i,k,0) = rhos(i,k,0)*rb(i,k)&
                        *sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)
                enddo
            enddo
        endif
    ENDDO

    elseif (vflag.eq.5) then
!
! TEST ONLY, vflag 5 and 6 are used to check the projection. Once we evaluate
! the new velocity fiels we evaluate the divergence to make sure that it is lower
! than the tolerance.
!
    do j=ys,ye
        do k=zs,ze
            do i=xs,xe
                uvel(i,k,j) = quarter*(q(i,k,j,1) + q(i-1,k,j,1)+&
                    q(i,k,j+1,1)+q(i-1,k,j+1,1))
                wvel(i,k,j) = quarter*(q(i,k,j,2) + q(i,k-1,j,2)+&
                    q(i,k,j+1,2)+q(i,k-1,j+1,2))
            enddo
        enddo
    enddo
    uvel(:,:,ye) = zero;wvel(:,:,ye) = zero

elseif (vflag.eq.6) then 
!
! TEST ONLY
!
    DO j=ys,ye
        do k=zs,ze
            do i=xs,xe
                rhogas = two*coe_rho/(f(i,k,j+1,1)+f(i,k,j,1))
                if(j == ny) rhogas = coe_rho/f(i,k,j,1)
                uu = uvel(i,k,j)
                ww = wvel(i,k,j)
                vv = q(i,k,j,3)
                vvel(i,k,j) = rhogas*(rb(i,k)*vec(i,k)+vv-uu*dphidx(i,k)-ww*dphidz(i,k))
            enddo
        enddo
        do k=zs-1,ze
            do i=xs-1,xe
                rhogas = two*coe_rho/(f(i,k,j,1)+f(i+1,k,j,1)) 
                uvel(i,k,j) = rhogas*q(i,k,j,1)
                rhogas = two*coe_rho/(f(i,k,j,1)+f(i,k+1,j,1)) 
                wvel(i,k,j) = rhogas*q(i,k,j,2)
            enddo
        enddo
    ENDDO

elseif (vflag.eq.7) then
!
!   VISUALIZATION ONLY, vflag ==7 is used for visualization purposes to put
!   the vectors back in the cell centered location so to visualize them using
!   either MATLAB OR ROCKETEER
!
    do k=zs,ze
        do i=xs,xe+1
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

!*********************************************************************
SUBROUTINE VEL_PRINT  !  [print rate]
  USE GLOBAL_DATA
  IMPLICIT NONE

!---------------------------------------------------------------------
! Local Variables
  INTEGER :: i, j
  REAL*8 :: uu, vv
!---------------------------------------------------------------------


  do i = 0, nx
     G_q(i,-1,3) = G_q(i,0,3)  
     do j = 0, ny
!
        uu = half*(G_q(i,j,1)+G_q(i-1,j,1))
        vv = half*(G_q(i,j,3)+G_q(i,j-1,3))
!        uu = G_q(i,j,1)
!        vv = G_q(i,j,3)

        G_rate(i,j,1) = uu
        G_rate(i,j,2) = vv

     end do
  end do

  do i = 0, nx
     do j = 0, ny

        G_q(i,j,1) = G_rate(i,j,1)
        G_q(i,j,2) = G_rate(i,j,2) 

     end do
  end do

  G_rate = 0.0d0



!----------------------------------------------------------     
  RETURN
END SUBROUTINE VEL_PRINT
!**************************************************************************
