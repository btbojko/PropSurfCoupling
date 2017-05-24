! ********************************************************************
SUBROUTINE EXP_RHS

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  include 'parallel.h'

!---------------------------------------------------------------------------------
! Local variables
  INTEGER ::  i, j, k, eqn, ys,ye,xs,xe,zs,ze, ip,kp,jp, col(4)
  REAL*8  :: dx1, dx2, dy1, dy2, dz1, dz2, dx1a, dx2a, dy1a, dy2a, dz1a, dz2a,term
  REAL*8  :: dxy1, dzy1,dlx1,dlx2,dlz1,dlz2,dly1,dly2,gx,gy,gxx,gyy,gxy,gz,gzz,gzy
  REAL*8  :: conv, exp_res, coe_exp
  REAL*8  :: diff, extra_term, pressure, poo
!---------------------------------------------------------------------------------

  poo = pressure(tcyc)

  ys = 1;ye=nyv(0)
  xs = drange(1);xe=nxv(0)
  zs = drange(3);ze=drange(4)

  dx1 = dx
  dx2 = dx*dx
  dy1 = dy
  dy2 = dy*dy
  dz1 = dz
  dz2 = dz*dz

  dx1a = 1.0d0/(2.0d0*dx1)
  dx2a = 1.0d0/(1.0d0*dx2)
  dy1a = 1.0d0/(2.0d0*dy1)
  dy2a = 1.0d0/(1.0d0*dy2)
  dz1a = 1.0d0/(2.0d0*dz1)
  dz2a = 1.0d0/(1.0d0*dz2)
  dxy1 = 1.0d0/(4.0d0*dx1*dy1)
  dzy1 = 1.0d0/(4.0d0*dz1*dy1)



!  dfdt(xs:xe,zs:ze,ys:ye,1:neqgas) = f(xs:xe,zs:ze,ys:ye,1:neqgas)

  do eqn =1,neqgas
  do i =xs,xe  
  do k =zs,ze  
  do j =ys,ye 
     if(eqn == 1) then
!        extra_term = zero
        extra_term = dpress_dt*dtx(i,k,j)/time_coeff/cp
     else
        extra_term = zero
     endif
     dqdt(i,k,j,eqn) = f(i,k,j,eqn) + extra_term
  enddo
  enddo
  enddo
  enddo

  if(neqmax == 1) RETURN

  eqn = neqmax
  coe_exp = merge(one-coe_dt, zero, coe_dt < one)
  do j = ys,ye
     do k = zs, ze
        do i = xs, xe

           dlx1 = 2.0d0*lambdas(i+1,k,j)*lambdas(i,k,j)/ &
                   (lambdas(i+1,k,j) + lambdas(i,k,j))
           dlx2 = 2.0d0*lambdas(i,k,j)*lambdas(i-1,k,j)/ &
                   (lambdas(i,k,j) + lambdas(i-1,k,j))
           dlz1 = 2.0d0*lambdas(i,k+1,j)*lambdas(i,k,j)/ &
                   (lambdas(i,k+1,j) + lambdas(i,k,j))
           dlz2 = 2.0d0*lambdas(i,k,j)*lambdas(i,k-1,j)/ &
                   (lambdas(i,k,j) + lambdas(i,k-1,j))
           dly1 = 2.0d0*lambdas(i,k,j+1)*lambdas(i,k,j)/ &
                   (lambdas(i,k,j+1) + lambdas(i,k,j))
           dly2 = 2.0d0*lambdas(i,k,j)*lambdas(i,k,j-1)/ &
                   (lambdas(i,k,j) + lambdas(i,k,j-1))

           gx = (f(i+1,k,j,eqn)-f(i-1,k,j,eqn))*dx1a
           gy = (f(i,k,j+1,eqn)-f(i,k,j-1,eqn))*dy1a
           gz = (f(i,k+1,j,eqn)-f(i,k-1,j,eqn))*dz1a

           

           gxx = (dlx1*(f(i+1,k,j,eqn)-f(i,k,j,eqn))&
                - dlx2*(f(i,k,j,eqn) - f(i-1,k,j,eqn)))*dx2a
           gyy = (dly1*(f(i,k,j+1,eqn) - f(i,k,j,eqn))&
                - dly2*(f(i,k,j,eqn) - f(i,k,j-1,eqn)))*dy2a
           gzz = (dlz1*(f(i,k+1,j,eqn)-f(i,k,j,eqn))&
                - dlz2*(f(i,k,j,eqn) - f(i,k-1,j,eqn)))*dz2a

           gxy = (+lambdas(i+1,k,j) * (f(i+1,k,j+1,eqn)&
                - f(i+1,k,j-1,eqn))&
                -lambdas(i-1,k,j) * (f(i-1,k,j+1,eqn)&
                - f(i-1,k,j-1,eqn))&
                +lambdas(i,k,j+1) * (f(i+1,k,j+1,eqn)&
                - f(i-1,k,j+1,eqn))&
                -lambdas(i,k,j-1) * (f(i+1,k,j-1,eqn)&
                - f(i-1,k,j-1,eqn)))*dxy1
           gzy = (+lambdas(i,k+1,j) * (f(i,k+1,j+1,eqn)&
                - f(i,k+1,j-1,eqn))&
                -lambdas(i,k-1,j) * (f(i,k-1,j+1,eqn)&
                - f(i,k-1,j-1,eqn))&
                +lambdas(i,k,j+1) * (f(i,k+1,j+1,eqn)&
                - f(i,k-1,j+1,eqn))&
                -lambdas(i,k,j-1) * (f(i,k+1,j-1,eqn)&
                - f(i,k-1,j-1,eqn)))*dzy1

           term = one+dphidx(i,k)**2+dphidz(i,k)**2

           exp_res = (vvel(i,k,0)*detady(j)*gy) &
                *rhos(i,k,j)/rhos(i,k,0)&
                + (gxx + gzz +&
                (term)*(&
                (detady(j)**2)*gyy+deta2dy2(j)*lambdas(i,k,j)*gy)&
                + dphidx(i,k)*detady(j)*gxy&
                + dphidz(i,k)*detady(j)*gzy&
                + (dphi2dx2(i,k)+dphi2dz2(i,k))*detady(j)*&
                lambdas(i,k,j)*gy)/cp
           dqdt(i,k,j,eqn) = coe_d * f(i,k,j,eqn) + coe_exp*dt*exp_res/rhos(i,k,j)
           rate(i,k,j,maxsize) = exp_res

        end do
     end do
  end do

!--------------------------------------------------------------------------
  RETURN
END SUBROUTINE EXP_RHS
!*********************************************************************


!**************************************************************************
SUBROUTINE EXP_RHSQ

  USE GLOBAL_DATA
  IMPLICIT NONE

!---------------------------------------------------------------------
! Local variables
  INTEGER ::  i, j, k, eqn
  INTEGER :: ys,ye,xs,xe,zs,ze, jp
  REAL*8  :: dx1, dx2, dy1, dy2, dz1, dz2, dx1a, dx2a, dy1a, dy2a, dz1a, dz2a,term
  REAL*8  :: ey,gx,gy,gz, conv, diff, gxx,gzz,gyy
  REAL*8  :: rhogas,dyy4,pressure,poo, prndtl_cp,coe_rho,coe_viscder
  REAL*8  :: firstcorr,corrflux,nu,nux,nuz,nuy
  REAL*8  :: s1,s2,s3,sig1,sig2,dya
  REAL*8  :: second_order,coe_exp
!---------------------------------------------------------------------
!
  prndtl_cp = pr/cp
  poo = PRESSURE(tcyc)
  coe_rho = dim_fact*poo

  dx1 = dx
  dy1 = dy
  dz1 = dz

  dx1a = 1.0d0/dx1
  dx2a = 1.0d0/(2.0d0*dx1)
  dy1a = 1.0d0/dy1
  dy2a = 1.0d0/(2.0d0*dy1)
  dz1a = 1.0d0/dz1
  dz2a = 1.0d0/(2.0d0*dz1)

  ys = drange(5)+1
  ye = nyv(0)
  zs = drange(3)
  ze = drange(4)
  xs = drange(1)
  xe = nxv(0)
  
  

  do j=0,nyv(0)+1
     do k=drange(3)-1,drange(4)+1
        do i=drange(1)-1,nxv(0)+1
           tmdpnt(i,k,j) = half*(f(i,k,j,1)+oldsoln(i,k,j,1))
           rhog(i,k,j) = coe_rho/tmdpnt(i,k,j)
        enddo
     enddo
  enddo
  

  oldsoln(:,:,:,1:ndim)=q(:,:,:,1:ndim)

  CALL LAMBDA(1)


  DO eqn = 1,ndim,2

     CALL VELOCITY(eqn)
     CALL VISC_DERIVATIVE(eqn)

     do j = ys, nyv(eqn)
        do k = zs, ze
           do i = xs, nxv(eqn)
              dfdt(i,k,j,eqn) =  q(i,k,j,eqn) &
                   + dtx(i,k,j)*dfdt(i,k,j,eqn) &
                   + dt*g_accel(eqn)  !gravity acceleration set in driver.f90 zero in most of cases
           end do
        end do
     end do

  ENDDO


!--- UPDATE the time-lagged pressure gradient and store it in the buffer rate
  do j=ys,nyv(1)
     dyy4 = quarter*detady(j)*dy1a
     do i=xs,nxv(1)
        do k=zs,ze
           rhogas = two*coe_rho/(f(i,k,j,1)+f(i+1,k,j,1))
           gy = (pold(i,k,j+1)+pold(i+1,k,j+1)-pold(i,k,j-1)-pold(i+1,k,j-1))*dyy4
           gx = (pold(i+1,k,j)-pold(i,k,j))*dx1a
           rate(i,k,j,1) = (dphidxa(i,k)*gy - gx)/rhogas
        enddo
     enddo
  enddo

  do j=ys,nyv(3)
     dya = detadya(j)*dy1a
     do i=xs,nxv(3)
        do k=zs,ze
           rhogas = two*coe_rho/(f(i,k,j+1,1)+f(i,k,j,1))
           gy = (pold(i,k,j+1)-pold(i,k,j))*dya
           rate(i,k,j,3) = - gy/rhogas
        enddo
     enddo
  enddo 

!---------------------------------------------------------------------
  RETURN
END SUBROUTINE EXP_RHSQ
!*********************************************************************


!*************************************************************************
SUBROUTINE VISCFLUX(eqn)  

  USE GLOBAL_DATA

  IMPLICIT NONE

!-------------------------------------------------------------------------  
  INTEGER,INTENT(IN) :: eqn
! Local variables
  INTEGER ::  i, j, k, ys,ye,xs,xe,zs,ze
  REAL*8  :: dx1, dx2, dy1, dy2, dz1, dz2, dx1a, dx4a, dy1a, dy4a, dz1a, dz4a,dyy,dyya
  REAL*8  :: gx,gy,gxx,gyy,gxy,gz,gzz,gzy,fxavg,fzavg
!-------------------------------------------------------------------------  

  ys = drange(5)-1
  ye = drange(6)+1
  if(ddrange(5).eq.0) ys = ys + 1
  if(ddrange(6).eq.ny) ye = ye - 1
  zs = drange(3)-1
  ze = drange(4)+1
  xs = drange(1)-1
  xe = drange(2)+1

! Initialize variables
  dx1 = dx
  dy1 = dy
  dz1 = dz
  
  dx1a = 1.0d0/dx1
  dx4a = 1.0d0/(4.0d0*dx1)
  dy1a = 1.0d0/dy1 
  dy4a = 1.0d0/(4.0d0*dy1)
  dz1a = 1.0d0/dz1 
  dz4a = 1.0d0/(4.0d0*dz1)
  
! THESE three loops could be (easily) recollected into one.
  if (eqn.eq.1) then
  do j = ys, ye
  dyy = dy4a*detady(j)
  dyya = dy1a*(detady(j)+detady(j+1))/2.0d0
  do k = zs, ze
  do i = xs, xe
              
!
! scalars @ i+1,j,k vectors i+1/2,j,k
!
     gx = (q(i+1,k,j,eqn)-q(i,k,j,eqn))*dx1a
     gy = (q(i,k,j+1,eqn)-q(i,k,j-1,eqn)+&
           q(i+1,k,j+1,eqn)-q(i+1,k,j-1,eqn))*dyy
     fdiff(i,k,j,1) = lambdag(i+1,k,j)*(gx-dphidx(i+1,k)*gy) 
!
! scalars @ i+1/2,j,k+1/2 vectors i,j,k+1/2
!
     gz = (q(i,k+1,j,eqn)-q(i,k,j,eqn))*dz1a
     gy = (q(i,k,j+1,eqn)-q(i,k,j-1,eqn)+&
           q(i,k+1,j+1,eqn)-q(i,k+1,j-1,eqn))*dyy
     fzavg = 0.5*(dphidza(i+1,k)+dphidza(i,k))
     fdiff(i,k,j,2) = lbdavg(i,k,j,1)*(gz-fzavg*gy)
!
! scalars @ i+1/2,j+1/2,k vectors i,j+1/2,k
!
     gx = (q(i+1,k,j,eqn)-q(i-1,k,j,eqn)+&
           q(i+1,k,j+1,eqn)-q(i-1,k,j+1,eqn))*dx4a
     gz = (q(i,k+1,j,eqn)-q(i,k-1,j,eqn)+&
           q(i,k+1,j+1,eqn)-q(i,k-1,j+1,eqn))*dz4a
     gy = (q(i,k,j+1,eqn)-q(i,k,j,eqn))*dyya
     fxavg = dphidxa(i,k)
     fzavg = 0.5*(dphidz(i+1,k)+dphidz(i,k))
     fdiff(i,k,j,3) = lbdavg(i,k,j,2)*(gy-fxavg*(gx-fxavg*gy)-&
                      fzavg*(gz-fzavg*gy))

   enddo
   enddo
   enddo

  elseif (eqn.eq.2) then

  do j = ys, ye
  dyy = dy4a*detady(j)
  dyya = dy1a*(detady(j)+detady(j+1))*half
  do k = zs, ze
  do i = xs, xe

     gx = (q(i+1,k,j,eqn)-q(i,k,j,eqn))*dx1a
     gy = (q(i,k,j+1,eqn)-q(i,k,j-1,eqn)+&
           q(i+1,k,j+1,eqn)-q(i+1,k,j-1,eqn))*dyy
     fxavg = 0.5*(dphidxa(i,k+1)+dphidxa(i,k))
     fdiff(i,k,j,1) = lbdavg(i,k,j,1)*(gx-fxavg*gy)
     gz = (q(i,k+1,j,eqn)-q(i,k,j,eqn))*dz1a
     gy = (q(i,k,j+1,eqn)-q(i,k,j-1,eqn)+&
           q(i,k+1,j+1,eqn)-q(i,k+1,j-1,eqn))*dyy
     fzavg = dphidz(i,k+1)
     fdiff(i,k,j,2) = lambdag(i,k+1,j)*(gz-fzavg*gy)
!    @   i,k+1/2,j+1/2
     gx = (q(i+1,k,j,eqn)-q(i-1,k,j,eqn)+&
           q(i+1,k,j+1,eqn)-q(i-1,k,j+1,eqn))*dx4a
     gz = (q(i,k+1,j,eqn)-q(i,k-1,j,eqn)+&
           q(i,k+1,j+1,eqn)-q(i,k-1,j+1,eqn))*dz4a
     gy = (q(i,k,j+1,eqn)-q(i,k,j,eqn))*dyya
     fxavg = 0.5*(dphidx(i,k+1)+dphidx(i,k))
     fzavg = dphidza(i,k)
     fdiff(i,k,j,3) = lbdavg(i,k,j,3)*(gy-fxavg*(gx-fxavg*gy)-&
                      fzavg*(gz-fzavg*gy))


   enddo
   enddo
   enddo

  elseif (eqn.eq.3) then
   
  do j = ys, ye
  dyy = dy4a*detadya(j)
  dyya = dy1a*detady(j+1)
  do k = zs, ze
  do i = xs, xe
     gx = (q(i+1,k,j,eqn)-q(i,k,j,eqn))*dx1a
     gy = (q(i,k,j+1,eqn)-q(i,k,j-1,eqn)+&
           q(i+1,k,j+1,eqn)-q(i+1,k,j-1,eqn))*dyy
     fxavg = dphidxa(i,k)
     fdiff(i,k,j,1) = lbdavg(i,k,j,2)*(gx-fxavg*gy)
     gz = (q(i,k+1,j,eqn)-q(i,k,j,eqn))*dz1a
     gy = (q(i,k,j+1,eqn)-q(i,k,j-1,eqn)+&
           q(i,k+1,j+1,eqn)-q(i,k+1,j-1,eqn))*dyy
     fzavg = dphidza(i,k)
     fdiff(i,k,j,2) = lbdavg(i,k,j,3)*(gz-fzavg*gy)
     gx = (q(i+1,k,j,eqn)-q(i-1,k,j,eqn)+&
           q(i+1,k,j+1,eqn)-q(i-1,k,j+1,eqn))*dx4a
     gz = (q(i,k+1,j,eqn)-q(i,k-1,j,eqn)+&
           q(i,k+1,j+1,eqn)-q(i,k-1,j+1,eqn))*dz4a
     gy = (q(i,k,j+1,eqn)-q(i,k,j,eqn))*dyya
     fxavg = dphidx(i,k)
     fzavg = dphidz(i,k)
     fdiff(i,k,j,3) = lambdag(i,k,j+1)*(gy-fxavg*(gx-fxavg*gy)-&
                      fzavg*(gz-fzavg*gy)) 

   enddo
   enddo
   enddo

   endif  !if (eqn.eq.1)
 
!--------------------------------------------------------------------
   RETURN
END SUBROUTINE VISCFLUX
!*********************************************************************

!*********************************************************************
SUBROUTINE VISC_DERIVATIVE(eqn)
!
! computes the terms due to the derivative of the viscosity
!
  USE GLOBAL_DATA
  IMPLICIT NONE

!--------------------------------------------------------------------
  ! Local Variables
  INTEGER :: i, j, k, ys,ye,xs,xe,zs,ze,eqn
  REAL *8 :: prndtl_cp,dx1a,dx2a,dx4a,dz1a,dz2a,dz4a,dy1a,dy2a,dy4a
  REAL *8 :: dx1,dz1,dy1,dyy,dy2,dy4,fxa,fza
  REAL *8 :: mx,mz,my,ux,uz,uy,wx,wy,vx,vz,vy,t1
!--------------------------------------------------------------------

  ! Initialize variables
  ys = drange(5)+1
  ye = nyv(eqn)
  zs = drange(3)
  ze = drange(4)
  xs = drange(1)
  xe = nxv(eqn)


!
! Initialize variables
!
  dx1 = dx
  dy1 = dy
  dz1 = dz
  
  dx1a = 1.0d0/dx1
  dx2a = 1.0d0/(2.0d0*dx1)
  dx4a = 1.0d0/(4.0d0*dx1)
  dy1a = 1.0d0/dy1 
  dy2a = 1.0d0/(2.0d0*dy1)
  dy4a = 1.0d0/(4.0d0*dy1)
  dz1a = 1.0d0/dz1 
  dz2a = 1.0d0/(2.0d0*dz1)
  dz4a = 1.0d0/(4.0d0*dz1)
  
  prndtl_cp = pr/cp
              
  if (eqn.eq.1) then

  do j = ys, ye
  dy4 = dy4a*detady(j)
  dy2 = dy2a*detady(j)
  do k = zs, ze
  do i = xs, xe

     my = prndtl_cp*(lambdag(i,k,j+1) - lambdag(i,k,j-1)+&
           lambdag(i+1,k,j+1)-lambdag(i+1,k,j-1))*dy4
     mx = prndtl_cp*(lambdag(i+1,k,j)-lambdag(i,k,j))*dx1a&
         - dphidxa(i,k)*my
     fza = 0.5*(dphidz(i,k)+dphidz(i+1,k))
     mz = prndtl_cp*(lambdag(i,k+1,j) - lambdag(i,k-1,j)+&
           lambdag(i+1,k+1,j)-lambdag(i+1,k-1,j))*dz4a&
         - fza*my

     vy = (q(i+1,k,j,3) + q(i,k,j,3) &
        -  q(i+1,k,j-1,3)-q(i,k,j-1,3))*dy2
     vx = (q(i+1,k,j,3) - q(i,k,j,3) &
        +  q(i+1,k,j-1,3)-q(i,k,j-1,3))*dx2a&
           - dphidxa(i,k)*vy

     wy = half*(q(i+1,k,j+1,2) + q(i+1,k-1,j+1,2) &
        - q(i+1,k,j-1,2) - q(i+1,k-1,j-1,2) & 
        + q(i,k,j+1,2)   + q(i,k-1,j+1,2)   &
        - q(i,k,j-1,2)   - q(i,k-1,j-1,2))*dy4
     wx = (q(i+1,k,j,2) - q(i,k,j,2) &
        +  q(i+1,k-1,j,2)-q(i,k-1,j,2))*dx2a&
           - dphidxa(i,k)*wy
     t1 = (q(i+1,k,j,2) + q(i,k,j,2)  &         !wz+vy
         - q(i+1,k-1,j,2)-q(i,k-1,j,2))*dz2a&
         - fza*wy + vy

     dfdt(i,k,j,eqn) = my*vx+mz*wx-mx*t1


  enddo
  enddo
  enddo

  elseif (eqn.eq.2) then

  do j = ys, ye
  dy4 = dy4a*detady(j)
  dy2 = dy2a*detady(j)
  do k = zs, ze
  do i = xs, xe

     fxa = 0.5*(dphidx(i,k)+dphidx(i,k+1))
     fza = dphidza(i,k)
     my = prndtl_cp*(lambdag(i,k,j+1) - lambdag(i,k,j-1)+&
           lambdag(i,k+1,j+1)-lambdag(i,k+1,j-1))*dy4
     mz = prndtl_cp*(lambdag(i,k+1,j)-lambdag(i,k,j))*dz1a&
         - fza*my
     mx = prndtl_cp*(lambdag(i+1,k,j) - lambdag(i-1,k,j)+&
           lambdag(i+1,k+1,j)-lambdag(i-1,k+1,j))*dx4a&
         - fxa*my

     vy = (q(i,k+1,j,3) + q(i,k,j,3) &
        -  q(i,k+1,j-1,3)-q(i,k,j-1,3))*dy2
     vz = (q(i,k+1,j,3) - q(i,k,j,3)+&
           q(i,k+1,j-1,3)-q(i,k,j-1,3))*dz2a&
         - fza*vy

     uy = half*(q(i-1,k,j+1,1) + q(i-1,k+1,j+1,1) &
              - q(i-1,k,j-1,1) - q(i-1,k+1,j-1,1) &
              + q(i,k,j+1,1)   + q(i,k+1,j+1,1)   &
              - q(i,k,j-1,1)   - q(i,k+1,j-1,1) )*dy4
!!     uy = (q(i,k,j+1,1) + q(i,k+1,j+1,1)-&
!!           q(i,k,j-1,1) - q(i,k+1,j-1,1))*dy4
     uz = (q(i,k+1,j,1) - q(i,k,j,1) &
         + q(i-1,k+1,j,1)-q(i-1,k,j,1))*dz2a&
         - fza*uy
     t1 = (q(i,k+1,j,1) + q(i,k,j,1)  &          !ux+vy
         - q(i-1,k+1,j,1)-q(i-1,k,j,1))*dx2a&
         - fxa*uy + vy

     dfdt(i,k,j,eqn) = mx*uz+my*vz-mz*t1

  enddo
  enddo
  enddo

  elseif (eqn.eq.3) then

  do j = ys, ye
  dy2 = dy2a*detady(j)
  dyy = dy1a*detadya(j)
  do k = zs, ze
  do i = xs, xe

     fxa = dphidx(i,k)
     fza = dphidz(i,k)
     my = prndtl_cp*(lambdag(i,k,j+1)-lambdag(i,k,j))*dyy
     mx = prndtl_cp*(lambdag(i+1,k,j+1) + lambdag(i+1,k,j)-&
             lambdag(i-1,k,j+1) - lambdag(i-1,k,j))*dx4a&
        - fxa*my
     mz = prndtl_cp*(lambdag(i,k+1,j+1) + lambdag(i,k+1,j)-&
             lambdag(i,k-1,j+1) - lambdag(i,k-1,j))*dz4a&
         - fza*my

     uy = (q(i,k,j+1,1) + q(i-1,k,j+1,1) &
         - q(i,k,j,1)  -  q(i-1,k,j,1))*dy2
     ux = (q(i,k,j+1,1) - q(i-1,k,j+1,1) &
         + q(i,k,j,1)  -  q(i-1,k,j,1))*dx2a&
         - fxa*uy
     wy = (q(i,k,j+1,2) + q(i,k-1,j+1,2) &
         - q(i,k,j,2)  -  q(i,k-1,j,2))*dy2
     t1 = (q(i,k,j+1,2) - q(i,k-1,j+1,2) &       !wz+ux
         + q(i,k,j,2)  -  q(i,k-1,j,2))*dz2a&
         - fza*wy + ux

     dfdt(i,k,j,eqn) = mx*uy+mz*wy-my*t1

  enddo
  enddo
  enddo

  endif

!---------------------------------------------------------------------
  RETURN
END SUBROUTINE VISC_DERIVATIVE
!*********************************************************************
