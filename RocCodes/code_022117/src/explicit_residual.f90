! ********************************************************************
SUBROUTINE EXP_RHS

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  include 'parallel.h'

!---------------------------------------------------------------------------------
! Local variables
  INTEGER ::  i, j, k, eqn,m,l, ys,ye,xs,xe,zs,ze, ip,kp,jp, col(4),uf(4),uDf(4)
  REAL*8  :: dx1, dx2, dy1, dy2, dz1, dz2, dx1a, dx2a, dy1a, dy2a, dz1a, dz2a,term
  REAL*8  :: dxy1, dzy1,dlmx1,dlmx2,dlmz1,dlmz2,dlmy1,dlmy2,gx,gy,gxx,gyy,gxy,gz,gzz,gzy
  REAL*8  :: conv,exp_res
  REAL*8  :: firstcorr,corrflux,s1,s2,s3,sig1,sig2,minmod,g2y
  REAL*8  :: diff,thetah_1,thetah_2,solid_dt
  REAL*8  :: gg(8),termxzsq,termxxzz,dlambdady,diffus
  REAL*8  :: t1,t31,t32,t33,DNM
!---------------------------------------------------------------------------------

  ys = 1
  ye = drange(6)
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)

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

  if(time_coeff == half) then
     thetah_1 = one
     thetah_2 = one
  else
     thetah_1 = one
     thetah_2 = zero
  endif
  coe_exp = thetah_1 + thetah_2

  mcomp = 1
  do i = drange(1),drange(2)
     do k = drange(3),drange(4)  

!independent of j
        termxzsq = dphidx(i,k)**2 + dphidz(i,k)**2
        termxxzz = dphi2dx2(i,k)  + dphi2dz2(i,k)

        do j = 1,drange(6)-1

           t1 = MQchi(j,mcomp)*termxzsq
           t31 = (one +  MQchipr(j,mcomp)*(MQphi(i,k)))
           t32 =  t31*t31
           t33 = t31*t32
           DNM = one/t33
!
           do eqn =1,neqgas
              diffus = lambdag(i,k,j)/lewis(eqn)
              dlambdady = dlgdy(i,k,j)*detady(j)/lewis(eqn)
              MQnonC(1:2) = diffus
              MQnonC(3) = diffus/t32*(one+t1*MQchi(j,mcomp))
              MQnonC(4) = -two*MQchi(j,mcomp)*dphidx(i,k)*diffus/t31
              MQnonC(5) = -two*MQchi(j,mcomp)*dphidz(i,k)*diffus/t31
              MQnonC(6) = dlgdx(i,k,j)/lewis(eqn) - MQchi(j,mcomp)*dphidx(i,k)*dlambdady/t31
              MQnonC(7) = dlgdz(i,k,j)/lewis(eqn) - MQchi(j,mcomp)*dphidz(i,k)*dlambdady/t31
              MQnonC(8) = DNM *( diffus*( t1* ( two*MQchipr(j,mcomp) + &
                   MQphi(i,k)*(two*MQchipr(j,mcomp)**2-MQchi(j,mcomp)*MQchisc(j,mcomp))) &
                   -MQchi(j,mcomp)*termxxzz*t32 - MQphi(i,k)*MQchisc(j,mcomp)) &
                   +dlambdady*t31*(one+MQchi(j,mcomp)*t1)&
                   -MQchi(j,mcomp)*t32*(dlgdx(i,k,j)*dphidx(i,k)+dlgdz(i,k,j)*dphidz(i,k))/lewis(eqn)   )

!change terms because the stretch mapping
              MQnonC(8) = MQnonC(8)*detady(j) + MQnonC(3)*deta2dy2(j)
              MQnonC(3) = MQnonC(3) *detady(j)**2
              MQnonC(4:5) = MQnonC(4:5)*detady(j)
!divide by the specific heat
              MQnonC = MQnonC/cp
!add convective terms
              MQnonC(6) = MQnonC(6) - Uvel(i,k,j)
              MQnonC(7) = MQnonC(7) - Wvel(i,k,j)
              MQnonC(8) = MQnonC(8) - Vvel(i,k,j)*detady(j)
!multiply by coe_exp*dt/rho
              MQnonC = MQnonC*dtx(i,k,j)*coe_exp


              gg(6) = (f(i+1,k,j,eqn)-f(i-1,k,j,eqn))*dx1a
              gg(7) = (f(i,k+1,j,eqn)-f(i,k-1,j,eqn))*dz1a
              gg(8) = (f(i,k,j+1,eqn)-f(i,k,j-1,eqn))*dy1a

              gg(1) = (f(i+1,k,j,eqn)-2.0d0*f(i,k,j,eqn)&
                   +f(i-1,k,j,eqn))*dx2a
              gg(2) = (f(i,k+1,j,eqn)-2.0d0*f(i,k,j,eqn)&
                   +f(i,k-1,j,eqn))*dz2a
              gg(3) = (f(i,k,j+1,eqn)-2.0d0*f(i,k,j,eqn)&
                   +f(i,k,j-1,eqn))*dy2a

              gg(4) = (f(i+1,k,j+1,eqn)-f(i+1,k,j-1,eqn)&
                   -f(i-1,k,j+1,eqn)+f(i-1,k,j-1,eqn))*dxy1
              gg(5) = (f(i,k+1,j+1,eqn)-f(i,k+1,j-1,eqn)&
                   -f(i,k-1,j+1,eqn)+f(i,k-1,j-1,eqn))*dzy1

              dfdt(i,k,j,eqn) = sum(gg*MQnonC) + coe_exp*rate(i,k,j,eqn)

!!>              do l = 1,ntrck
!!>                 if(i-itrck(l) == 0 .and. j - jtrck(l) == 0.and. k - ktrck(l) == 0.and. myid - mtrck(l) == 0) then
!!>                    print*,eqn,myid,i,k,j,dfdt(i,k,j,eqn), sum(gg*MQnonC),'|||||',gg(1:8),'|||||',MQnonC(1:8),dtx(i,k,j)
!!>                 end if
!!>              end do

           enddo
        enddo
     enddo
  enddo

!
  if(inewton == 1 .and. NNewton > 1) then
     rhsOLD(:,:,:,1:neqgas) = thetah_2*dfdt(:,:,:,1:neqgas)
  elseif(inewton > 1) then
     dfdt(:,:,:,1:neqgas) = thetah_1*dfdt(:,:,:,1:neqgas)+ rhsOLD(:,:,:,1:neqgas)
  endif


  if(ipack<0 .and. mod(ncyc,writemod) == 2) then
     write(*,*) ncyc,'MAX res',maxval(abs(dfdt(drange(1):drange(2), drange(3):drange(4) ,1:drange(6)-1,1:neqgas))),&
          &maxloc(abs(dfdt(drange(1):drange(2), drange(3):drange(4) ,1:drange(6)-1,1:neqgas)))+&
          &(/drange(1), drange(3) ,1,1/) -1
  end if

  uf = Ubound(f)
  uf(4) = neqgas
  uDf = Ubound(dfdt)
  uDf(4) = neqgas

  call CORR_FLUX(f,lbound(f),uf,dfdt,lbound(dfdt),uDf)

  eqn = neqmax
  if(coe_dt == half) then
     thetah_1 = coe_dt
     thetah_2 = coe_dt
  else
     thetah_1 = coe_dt
     thetah_2 = zero
  endif

  coe_exp = thetah_1 + thetah_2
  solid_dt = coe_exp*dt

  do j = 1,drange(6)-1
     do k = drange(3),drange(4)  
        do i = drange(1),drange(2)

           dlmx1 = 2.0d0*lambdas(i+1,k,j)*lambdas(i,k,j)/ &
                (lambdas(i+1,k,j) + lambdas(i,k,j))
           dlmx2 = 2.0d0*lambdas(i,k,j)*lambdas(i-1,k,j)/ &
                (lambdas(i,k,j) + lambdas(i-1,k,j))
           dlmz1 = 2.0d0*lambdas(i,k+1,j)*lambdas(i,k,j)/ &
                (lambdas(i,k+1,j) + lambdas(i,k,j))
           dlmz2 = 2.0d0*lambdas(i,k,j)*lambdas(i,k-1,j)/ &
                (lambdas(i,k,j) + lambdas(i,k-1,j))
           dlmy1 = 2.0d0*lambdas(i,k,j+1)*lambdas(i,k,j)/ &
                (lambdas(i,k,j+1) + lambdas(i,k,j))
           dlmy2 = 2.0d0*lambdas(i,k,j)*lambdas(i,k,j-1)/ &
                (lambdas(i,k,j) + lambdas(i,k,j-1))

           gx = (f(i+1,k,j,eqn)-f(i-1,k,j,eqn))*dx1a
           gy = (f(i,k,j+1,eqn)-f(i,k,j-1,eqn))*dy1a
           gz = (f(i,k+1,j,eqn)-f(i,k-1,j,eqn))*dz1a

           gxx = (dlmx1*(f(i+1,k,j,eqn)-f(i,k,j,eqn))&
                - dlmx2*(f(i,k,j,eqn) - f(i-1,k,j,eqn)))*dx2a
           gyy = (dlmy1*(f(i,k,j+1,eqn) - f(i,k,j,eqn))&
                - dlmy2*(f(i,k,j,eqn) - f(i,k,j-1,eqn)))*dy2a
           gzz = (dlmz1*(f(i,k+1,j,eqn)-f(i,k,j,eqn))&
                - dlmz2*(f(i,k,j,eqn) - f(i,k-1,j,eqn)))*dz2a

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

           dfdt(i,k,j,eqn) = solid_dt*exp_res/rhos(i,k,j)

        end do
     end do
  end do
  if(inewton == 1 .and. NNewton > 1) then
     rhsOLD(:,:,:,neqmax) = thetah_2*dfdt(:,:,:,neqmax)
  elseif(inewton > 1) then
     dfdt(:,:,:,neqmax) = thetah_1*dfdt(:,:,:,neqmax)+ rhsOLD(:,:,:,neqmax)
  endif

!--------------------------------------------------------------------------
  RETURN
END SUBROUTINE EXP_RHS
!*********************************************************************

!****************************************************************************************
SUBROUTINE CORR_FLUX(Fvec,lf,uf,DFvec,lDf,uDf)

  USE GLOBAL_DATA

  IMPLICIT NONE

  include 'parallel.h'

!----------------------------------------------------------------------------------------
! Local variables
  INTEGER :: i,j,k, m, ys,ye,xs,xe,zs,ze
  INTEGER :: neqn, ip,kp,jp, idir
  INTEGER :: lf(4),uf(4),lDf(4),uDf(4)
  INTEGER,DIMENSION(3) :: ikj,ikj_P, ikj_P1,ikj_M1,ikj_M2,ikjP,ikjM
  REAL*8  :: corrflux,s1,s2,s3,sig1,sig2,nuy
  REAL*8  :: TERM_FD,invcp,dtcp
  REAL*8  :: Fvec(lf(1):uf(1),lf(2):uf(2),lf(3):uf(3),lf(4):uf(4))
  REAL*8  :: DFvec(lDf(1):uDf(1),lDf(2):uDf(2),lDf(3):uDf(3),lDf(4):uDf(4))
!----------------------------------------------------------------------------------------

!
!  to improve performance in this 3-D version use TVD only in y direction
!
  neqn = uf(4)

  ys = 1
  ye = drange(6)-1
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)


  DO j = ys, ye
     do k = zs, ze
        do i = xs, xe

           ikj = (/i,k,j/)
           do idir = 1,ndim
 
              nuy = dconv(i,k,j,idir)

              if(nuy < 1.d-9) CYCLE

              ikj_P = ikj;
              ikj_P(idir) = ikj(idir)+ iconv(i,k,j,idir)
              ikj_P1 = ikj_P; ikj_P1(idir) = ikj_P(idir)+1
              ikj_M1 = ikj_P; ikj_M1(idir) = ikj_P(idir)-1
              ikj_M2 = ikj_P; ikj_M2(idir) = ikj_P(idir)-2

              ikjP = ikj;ikjP(idir) = ikj(idir)+1
              ikjM = ikj;ikjM(idir) = ikj(idir)-1

              do m = 1,neqn
                 TERM_FD = (Fvec(ikjP(1),ikjP(2),ikjP(3),m) - two*Fvec(ikj(1),ikj(2),ikj(3),m) &
                      + Fvec(ikjM(1),ikjM(2),ikjM(3),m))

                 if(.not. skipTVD) then
                    s1 = (Fvec(ikj_P1(1),ikj_P1(2),ikj_P1(3),m) - Fvec(ikj_P(1),ikj_P(2),ikj_P(3),m))
                    s2 = (Fvec(ikj_P(1),ikj_P(2),ikj_P(3),m) - Fvec(ikj_M1(1),ikj_M1(2),ikj_M1(3),m))
                    s3 = (Fvec(ikj_M1(1),ikj_M1(2),ikj_M1(3),m) - Fvec(ikj_M2(1),ikj_M2(2),ikj_M2(3),m))

                    sig1 = (sign(half,s1)+sign(half,s2))*min(abs(s1),abs(s2))
                    sig2 = (sign(half,s2)+sign(half,s3))*min(abs(s2),abs(s3))

                    TERM_FD = TERM_FD - (sig1-sig2)
                 endif

                 corrflux =  nuy * TERM_FD

                 DFvec(i,k,j,m) = DFvec(i,k,j,m) + corrflux*coe_exp
             
              enddo

           enddo

        enddo
     enddo
  ENDDO
!------------------------------------------------------------------------------------
  RETURN 
END SUBROUTINE CORR_FLUX
!****************************************************************************************


!**************************************************************************
SUBROUTINE EXP_RHSQ(eqn)

  USE GLOBAL_DATA
  IMPLICIT NONE

  !---------------------------------------------------------------------
  INTEGER, INTENT(IN) :: eqn
  ! Local variables
  INTEGER ::  i, j, k, ys,ye,xs,xe,zs,ze,lf(4),uf(4),lDf(4),uDf(4)
  REAL*8  :: dx1, dx2, dy1, dy2, dz1, dz2, dx1a, dx2a, dy1a, dy2a, dz1a, dz2a
  REAL*8  :: ey,dyy,gx,gy,gz, conv, diff, gyy
  REAL*8  :: rhogas,dya4,pressure,poo, prndtl_cp
  REAL*8  :: coe_rho,coe_viscder
  !---------------------------------------------------------------------
  !
  prndtl_cp = pr/cp
  poo = PRESSURE(tcyc(2))
  coe_rho = two*dim_fact*poo

  dx1 = dx
  dy1 = dy
  dz1 = dz


  !-- dx1a first derivative, dx2a second derivative
  dx1a = 1.0d0/(2.0d0*dx1)
  dx2a = 1.0d0/(1.0d0*dx1)
  dy1a = 1.0d0/(2.0d0*dy1)
  dy2a = 1.0d0/(1.0d0*dy1)
  dz1a = 1.0d0/(2.0d0*dz1)
  dz2a = 1.0d0/(1.0d0*dz1)

  ys = drange(5)+1
  ye = drange(6)-1
  zs = drange(3)
  ze = drange(4)
  xs = drange(1)
  xe = drange(2)



  CALL VISC_DERIVATIVE(eqn)
  CALL VISCFLUX(eqn)


  coe_viscder = merge(one/time_coeff, zero, ncyc > ncyc_init .and. .not. skipViscDeriv)
  coe_exp = merge(one, zero, time_coeff == half)

  do j = ys, ye
     ey = merge(detady(j),detadya(j),eqn <3)
     dyy=dy2a*ey
     do k = zs, ze
        do i = xs, xe
           if(eqn == 1)MQterm(i,k,1) = (one +  MQchipr(j,1)*(MQphi(i+1,k)+MQphi(i,k))*half)
           if(eqn == 2)MQterm(i,k,2) = (one +  MQchipr(j,1)*(MQphi(i,k+1)+MQphi(i,k))*half)
           if(eqn == 3)MQterm(i,k,3) = (one +  MQchipr(j,2)*MQphi(i,k))

           gx = (q(i+1,k,j,eqn)-q(i-1,k,j,eqn))*dx1a
           gz = (q(i,k+1,j,eqn)-q(i,k-1,j,eqn))*dz1a
           gy = (q(i,k,j+1,eqn)-q(i,k,j-1,eqn))*dy1a

           conv = - vvel(i,k,j)*gy*ey &
                - uvel(i,k,j)*gx - wvel(i,k,j)*gz

           diff = (dff(i,k,j,1)-dff(i-1,k,j,1))*dx2a&
                +(dff(i,k,j,2)-dff(i,k-1,j,2))*dz2a&
                +(dff(i,k,j,3)-dff(i,k,j-1,3))*dyy

           dfdt(i,k,j,eqn) =  dtx(i,k,j)*dfdt(i,k,j,eqn)*coe_viscder    &
                &          +  q(i,k,j,eqn)                              &
                &          +  coe_exp*dtx(i,k,j)*(conv + Prndtl_cp*diff/MQterm(i,k,eqn))
        end do
     end do
  end do

!---
  IF(time_coeff == half) THEN

     select case (eqn)

     case(1) 
        do j=ys,ye
           dya4 = quarter*detady(j)*dy1a
           do k = zs-1, ze+1
              do i = xs-1, xe+1
                 MQterm(i,k,0) = (one +  MQchipr(j,1)*(MQphi(i,k)))
                 MQterm(i,k,1) = (one +  MQchipr(j,1)*(MQphi(i+1,k)+MQphi(i,k))*half)
              enddo
           enddo
           do i=xs,xe
              do k=zs,ze
                 if(SkipVarDensity) then
                    rhogas = rhoRef
                 else
                    rhogas = coe_rho/(f(i,k,j,1)+f(i+1,k,j,1))
                 endif
                 gy = (MQchi(j+1,1)*(pold(i,k,j+1)+pold(i+1,k,j+1)) - &
                      &MQchi(j-1,1)*(pold(i,k,j-1)+pold(i+1,k,j-1)))*dya4
                 gx = (MQterm(i+1,k,0)*pold(i+1,k,j)-MQterm(i,k,0)*pold(i,k,j))*dx1a

                 rate(i,k,j,1) = (dphidxa(i,k)*gy - gx)/rhogas/MQterm(i,k,1)
              enddo
           enddo
        enddo

     case(2)  
        do j=ys,ye
           dya4 = quarter*detady(j)*dy1a
           do k = zs-1, ze+1
              do i = xs-1, xe+1
                 MQterm(i,k,0) = (one +  MQchipr(j,1)*(MQphi(i,k)))
                 MQterm(i,k,2) = (one +  MQchipr(j,1)*(MQphi(i,k+1)+MQphi(i,k))*half)
              enddo
           enddo
           do i=xs,xe
              do k=zs,ze

                 if(SkipVarDensity) then
                    rhogas = rhoRef
                 else
                    rhogas = coe_rho/(f(i,k,j,1)+f(i,k+1,j,1))
                 endif

                 gy = (MQchi(j+1,1)*(pold(i,k,j+1)+pold(i,k+1,j+1))-&
                      &MQchi(j+1,1)*(pold(i,k,j-1)+pold(i,k+1,j-1)))*dya4
                 gz = (MQterm(i,k+1,0)*pold(i,k+1,j)-MQterm(i,k,0)*pold(i,k,j))*dz1a

                 rate(i,k,j,2) = (dphidza(i,k)*gy - gz)/rhogas/MQterm(i,k,2)
              enddo
           enddo
        enddo


     case(3)  
        do j=ys,ye
           dyy = detadya(j)*dy1a
           do k = zs-1, ze+1
              do i = xs-1, xe+1
                 MQterm(i,k,0) = (one +  MQchipr(j,1)*(MQphi(i,k)))
                 MQterm(i,k,3) = (one +  MQchipr(j,2)*(MQphi(i,k)))
              enddo
           enddo
           do i=xs,xe
              do k=zs,ze
                 if(SkipVarDensity) then
                    rhogas = rhoRef
                 else
                    rhogas = coe_rho/(f(i,k,j+1,1)+f(i,k,j,1))
                 endif
                 gy = (pold(i,k,j+1)-pold(i,k,j))*dyy
                 rate(i,k,j,3) = - gy/rhogas/MQterm(i,k,3)
              enddo
           enddo
        enddo

     END select

     dqdt(:,:,:,eqn) = dfdt(:,:,:,eqn) +  rate(:,:,:,eqn)
!!!---Initialize the pressure perturbation
     p = zero

  ELSE

     pold = zero
     dqdt(:,:,:,eqn) = dfdt(:,:,:,eqn)

  ENDIF

  if(time_coeff == one) RETURN


  uf = Ubound(q)
  uf(4) = 1
  lf = lbound(q)
  lf(4) = 1
  uDf = Ubound(dqdt)
  uDf(4) = 1
  lDf = Lbound(dqdt)
  lDf(4) = 1

  call CORR_FLUX(q(lf(1),lf(2),lf(3),eqn),lf,uf,dqdt(lDf(1),lDf(2),lDf(3),eqn),lDf,uDf)


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
  REAL*8  :: dx1, dx2, dy1, dy2, dz1, dz2, dx1a, dx4a, dy1a, dy4a, dz1a, dz4a,dyy,dyyh
  REAL*8  :: gx,gy,gxx,gyy,gxy,gz,gzz,gzy,fxavg,fzavg
!-------------------------------------------------------------------------  

  ys = drange(5)-1
  ye = drange(6)+1
  if(dSerialRange(5).eq.0) ys = ys + 1
  if(dSerialRange(6).eq.ny) ye = ye - 1
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

  select case (eqn)

  case(1)

     mcomp  = 1;
     mgradj = 2
     do j = ys, ye
        dyy = dy4a*detady(j)
        dyyh = dy1a*detadya(j)

        do k = zs-1, ze+1
           do i = xs-1, xe+1
              do Mindx = 1,2
                 MQterm(i,k,Mindx) = (one +  MQchipr(j,Mindx)*(MQphi(i+1,k)+MQphi(i,k))*half)
              enddo
           enddo
        enddo

        do k = zs, ze
           do i = xs, xe

!
! scalars @ i+1,j,k vectors i+1/2,j,k
!
              gx = (MQterm(i+1,k,mcomp)*q(i+1,k,j,eqn)-MQterm(i,k,mcomp)*q(i,k,j,eqn))*dx1a
              gy = ( (q(i,k,j+1,eqn)+q(i+1,k,j+1,eqn)) * MQchi(j+1,mcomp)&
                   - (q(i,k,j-1,eqn)+q(i+1,k,j-1,eqn)) * MQchi(j-1,mcomp))*dyy
              dff(i,k,j,1) = lambdag(i+1,k,j)*(gx-dphidx(i+1,k)*gy) 
!
! scalars @ i+1/2,j,k+1/2 vectors i,j,k+1/2
!
              gz = (MQterm(i,k+1,mcomp)*q(i,k+1,j,eqn)-MQterm(i,k,mcomp)*q(i,k,j,eqn))*dz1a
              gy = ((q(i,k,j+1,eqn)+q(i,k+1,j+1,eqn)) * MQchi(j+1,mcomp)&
                   -(q(i,k,j-1,eqn)+q(i,k+1,j-1,eqn)) * MQchi(j-1,mcomp))*dyy
              fzavg = 0.5*(dphidza(i+1,k)+dphidza(i,k))
              dff(i,k,j,2) = lbdavg(i,k,j,1)*(gz-fzavg*gy)
!
! scalars @ i+1/2,j+1/2,k vectors i,j+1/2,k
!
              gx   = MQchi(j,mgradj)*((q(i+1,k,j,eqn)+q(i+1,k,j+1,eqn))*MQterm(i+1,k,mgradj)&
                   -(q(i-1,k,j,eqn)+q(i-1,k,j+1,eqn))*MQterm(i-1,k,mgradj))*dx4a
              gz   = MQchi(j,mgradj)*((q(i,k+1,j,eqn)+q(i,k+1,j+1,eqn))*MQterm(i,k+1,mgradj)&
                   -(q(i,k-1,j,eqn)+q(i,k-1,j+1,eqn))*MQterm(i,k-1,mgradj))*dz4a
              MQgy = MQchi(j,mgradj)*(MQchi(j+1,mcomp)*q(i,k,j+1,eqn)-MQchi(j,mcomp)*q(i,k,j,eqn))*dyyh
              gy = (q(i,k,j+1,eqn)-q(i,k,j,eqn))*dyyh
              fxavg = dphidxa(i,k)
              fzavg = 0.5*(dphidz(i+1,k)+dphidz(i,k))
              dff(i,k,j,3) = lbdavg(i,k,j,2)*(gy+MQgy*(fxavg**2+fzavg**2)-fxavg*gx-&
                   &fzavg*gz)/MQterm(i,k,mgradj)

           enddo
        enddo
     enddo

  case(2)

     mcomp = 1;
     mgradj = 2
     do j = ys, ye
        dyy = dy4a*detady(j)
        dyyh = dy1a*detadya(j)

        do k = zs-1, ze+1
           do i = xs-1, xe+1
              do Mindx = 1,2
                 MQterm(i,k,Mindx) = (one +  MQchipr(j,Mindx)*(MQphi(i,k+1)+MQphi(i,k))*half)
              enddo
           enddo
        enddo

        do k = zs, ze
           do i = xs, xe

              gx = (MQterm(i+1,k,mcomp)*q(i+1,k,j,eqn)-MQterm(i,k,mcomp)*q(i,k,j,eqn))*dx1a
              gy = ( (q(i,k,j+1,eqn)+q(i+1,k,j+1,eqn)) * MQchi(j+1,mcomp)&
                   - (q(i,k,j-1,eqn)+q(i+1,k,j-1,eqn)) * MQchi(j-1,mcomp))*dyy
              fxavg = 0.5*(dphidxa(i,k+1)+dphidxa(i,k))
              dff(i,k,j,1) = lbdavg(i,k,j,1)*(gx-fxavg*gy)

              gz = (MQterm(i,k+1,mcomp)*q(i,k+1,j,eqn)-MQterm(i,k,mcomp)*q(i,k,j,eqn))*dz1a
              gy = ((q(i,k,j+1,eqn)+q(i,k+1,j+1,eqn)) * MQchi(j+1,mcomp)&
                   -(q(i,k,j-1,eqn)+q(i,k+1,j-1,eqn)) * MQchi(j-1,mcomp))*dyy
              fzavg = dphidz(i,k+1)
              dff(i,k,j,2) = lambdag(i,k+1,j)*(gz-fzavg*gy)
!    @   i,k+1/2,j+1/2
!
              gx   = MQchi(j,mgradj)*((q(i+1,k,j,eqn)+q(i+1,k,j+1,eqn))*MQterm(i+1,k,mgradj)&
                   -(q(i-1,k,j,eqn)+q(i-1,k,j+1,eqn))*MQterm(i-1,k,mgradj))*dx4a
              gz   = MQchi(j,mgradj)*((q(i,k+1,j,eqn)+q(i,k+1,j+1,eqn))*MQterm(i,k+1,mgradj)&
                   -(q(i,k-1,j,eqn)+q(i,k-1,j+1,eqn))*MQterm(i,k-1,mgradj))*dz4a
              MQgy = MQchi(j,mgradj)*(MQchi(j+1,mcomp)*q(i,k,j+1,eqn)-MQchi(j,mcomp)*q(i,k,j,eqn))*dyyh
              gy = (q(i,k,j+1,eqn)-q(i,k,j,eqn))*dyyh

              fxavg = 0.5*(dphidx(i,k+1)+dphidx(i,k))
              fzavg = dphidza(i,k)
              dff(i,k,j,3) = lbdavg(i,k,j,3)*(gy+MQgy*(fxavg**2+fzavg**2)-fxavg*gx-&
                   &fzavg*gz)/MQterm(i,k,mgradj)
           enddo
        enddo
     enddo

  case(3)

     mcomp = 2;
     mgradj = 1
     do j = ys, ye
        dyy = dy4a*detadya(j)
        dyyh = dy1a*detady(j+1)

        do k = zs-1, ze+1
           do i = xs-1, xe+1
              MQterm(i,k,2) = (one +  MQchipr(j,2)*MQphi(i,k))
              MQterm(i,k,1) = (one +  MQchipr(j+1,1)*MQphi(i,k))
           enddo
        enddo

        do k = zs, ze
           do i = xs, xe
              gx = (MQterm(i+1,k,mcomp)*q(i+1,k,j,eqn)-MQterm(i,k,mcomp)*q(i,k,j,eqn))*dx1a
              gy = ( (q(i,k,j+1,eqn)+q(i+1,k,j+1,eqn)) * MQchi(j+1,mcomp)&
                   - (q(i,k,j-1,eqn)+q(i+1,k,j-1,eqn)) * MQchi(j-1,mcomp))*dyy
              fxavg = dphidxa(i,k)
              dff(i,k,j,1) = lbdavg(i,k,j,2)*(gx-fxavg*gy)

              gz = (MQterm(i,k+1,mcomp)*q(i,k+1,j,eqn)-MQterm(i,k,mcomp)*q(i,k,j,eqn))*dz1a
              gy = ((q(i,k,j+1,eqn)+q(i,k+1,j+1,eqn)) * MQchi(j+1,mcomp)&
                   -(q(i,k,j-1,eqn)+q(i,k+1,j-1,eqn)) * MQchi(j-1,mcomp))*dyy
              fzavg = dphidza(i,k)
              dff(i,k,j,2) = lbdavg(i,k,j,3)*(gz-fzavg*gy)
!
              gx   = MQchi(j+1,mgradj)*((q(i+1,k,j,eqn)+q(i+1,k,j+1,eqn))*MQterm(i+1,k,mgradj)&
                   -(q(i-1,k,j,eqn)+q(i-1,k,j+1,eqn))*MQterm(i-1,k,mgradj))*dx4a
              gz   = MQchi(j+1,mgradj)*((q(i,k+1,j,eqn)+q(i,k+1,j+1,eqn))*MQterm(i,k+1,mgradj)&
                   -(q(i,k-1,j,eqn)+q(i,k-1,j+1,eqn))*MQterm(i,k-1,mgradj))*dz4a
              MQgy = MQchi(j+1,mgradj)*(MQchi(j+1,mcomp)*q(i,k,j+1,eqn)-MQchi(j,mcomp)*q(i,k,j,eqn))*dyyh
              gy = (q(i,k,j+1,eqn)-q(i,k,j,eqn))*dyyh
              fxavg = dphidx(i,k)
              fzavg = dphidz(i,k)
              dff(i,k,j,3) = lambdag(i,k,j+1)*(gy+MQgy*(fxavg**2+fzavg**2)-fxavg*gx-&
                   &fzavg*gz)/MQterm(i,k,mgradj)

           enddo
        enddo
     enddo

  end select

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
  ys = drange(5)
  ye = drange(6)
  if(dSerialRange(5).eq.0) ys = ys + 1
  if(dSerialRange(6).eq.ny) ye = ye - 1
  zs = drange(3)
  ze = drange(4)
  xs = drange(1)
  xe = drange(2)


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
  
  mcomp = merge(1,2,eqn<=2)
              
  if (eqn.eq.1) then

  do j = ys, ye
  do k = zs, ze
  do i = xs, xe

     dy4 = dy4a*detady(j)/(one+MQchipr(j,mcomp)*MQphi(i,k))
     dy2 = dy2a*detady(j)/(one+MQchipr(j,mcomp)*MQphi(i,k))
     fza = half*(dphidz(i,k)+dphidz(i+1,k))*MQchi(j,mcomp)
     fxa = dphidxa(i,k)*MQchi(j,mcomp)

     my = prndtl_cp*(lambdag(i,k,j+1) - lambdag(i,k,j-1)+&
           lambdag(i+1,k,j+1)-lambdag(i+1,k,j-1))*dy4
     mx = prndtl_cp*(lambdag(i+1,k,j)-lambdag(i,k,j))*dx1a&
         - fxa*my
     mz = prndtl_cp*(lambdag(i,k+1,j) - lambdag(i,k-1,j)+&
           lambdag(i+1,k+1,j)-lambdag(i+1,k-1,j))*dz4a&
         - fza*my

     vy = (q(i+1,k,j,3) + q(i,k,j,3) &
        -  q(i+1,k,j-1,3)-q(i,k,j-1,3))*dy2
     vx = (q(i+1,k,j,3) - q(i,k,j,3) &
        +  q(i+1,k,j-1,3)-q(i,k,j-1,3))*dx2a&
           - fxa*vy

     wy = half*(q(i+1,k,j+1,2) + q(i+1,k-1,j+1,2) &
        - q(i+1,k,j-1,2) - q(i+1,k-1,j-1,2) & 
        + q(i,k,j+1,2)   + q(i,k-1,j+1,2)   &
        - q(i,k,j-1,2)   - q(i,k-1,j-1,2))*dy4
     wx = (q(i+1,k,j,2) - q(i,k,j,2) &
        +  q(i+1,k-1,j,2)-q(i,k-1,j,2))*dx2a&
           - fxa*wy
     t1 = (q(i+1,k,j,2) + q(i,k,j,2)  &         !wz+vy
         - q(i+1,k-1,j,2)-q(i,k-1,j,2))*dz2a&
         - fza*wy + vy

     dfdt(i,k,j,eqn) = my*vx+mz*wx-mx*t1
     
  enddo
  enddo
  enddo

  elseif (eqn.eq.2) then

  do j = ys, ye
  do k = zs, ze
  do i = xs, xe

     dy4 = dy4a*detady(j)/(one+MQchipr(j,mcomp)*MQphi(i,k))
     dy2 = dy2a*detady(j)/(one+MQchipr(j,mcomp)*MQphi(i,k))

     fxa = half*(dphidx(i,k)+dphidx(i,k+1))*MQchi(j,mcomp)
     fza = dphidza(i,k)*MQchi(j,mcomp)

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
  do k = zs, ze
  do i = xs, xe

     
     dy2 = dy2a*detadya(j)/(one+MQchipr(j,mcomp)*MQphi(i,k))
     dyy = dy1a*detadya(j)/(one+MQchipr(j,mcomp)*MQphi(i,k))


     fxa = dphidx(i,k)*MQchi(j,mcomp)
     fza = dphidz(i,k)*MQchi(j,mcomp)

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
! *******************************************************************
SUBROUTINE CONTINUITYSOURCE

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  include 'parallel.h'

!--------------------------------------------------------------------------------- 
! Local variables 
  INTEGER ::  i, j, k, eqn
  REAL*8  :: pressure
!--------------------------------------------------------------------------------- 

  rhogasCont = cshift(rhogasCont,-1,4)
  dtCont = cshift(dtCont,-1)
  
  dtCont(1) = timestep
  do k=drange(3),drange(4)
     do i=drange(1),drange(2)
        do j=drange(5),drange(6)
           rhogasCont(i,k,j,1) = dim_fact*PRESSURE(tcyc(2))/f(i,k,j,1)
        enddo
     enddo
  enddo

  rhogasCont(:,:,jRFLO(1):ubound(rhogasCont,3),:) = zero

  return
end SUBROUTINE CONTINUITYSOURCE
! ************************************************************** 


! **************************************************************
SUBROUTINE CONTINUITYSOURCEINIT

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  include 'parallel.h'

!---------------------------------------------------------------------------------
! Local variables  
  INTEGER ::  i, j, k, eqn
  REAL*8  :: pressure
!---------------------------------------------------------------------------------

  if (irestart == 0) then
     dtCont = (/timeStep,timeStep /)
  else
     dtCont = (/timeStep,timeStep  /)
  endif
  do k=drange(3),drange(4)
     do i=drange(1),drange(2)
        do j=drange(5),drange(6)
           rhogasCont(i,k,j,1:3) = dim_fact*PRESSURE(tcyc(2))/f(i,k,j,1)
        enddo
     enddo
  enddo
  
  return
end SUBROUTINE CONTINUITYSOURCEINIT


! **************************************************************                                  
SUBROUTINE MidPointDensityConductivity

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  include 'parallel.h'

  !---------------------------------------------------------------------------------                
  ! Local variables                                                                                 
  INTEGER ::  i, j, k, eqn
  REAL*8  :: pressure,coe_rho
  !---------------------------------------------------------------------------------                
  if(.not. doprojection) return

  coe_rho = dim_fact*PRESSURE(tcyc(2))

  coe_rho = dim_fact*press

  do j=0,ny
     do k=drange(3)-1,drange(4)+1
        do i=drange(1)-1,drange(2)+1

           tmdpnt(i,k,j) = half*(f(i,k,j,1)+oldsoln(i,k,j,1))
           rhog(i,k,j) = coe_rho/tmdpnt(i,k,j)

        enddo
     enddo
  enddo
  oldsoln(:,:,:,maxsize+1:maxsize+ndim)=q(:,:,:,1:ndim)
  CALL LAMBDA(1)

end SUBROUTINE MidPointDensityConductivity

! **************************************************************                                  
SUBROUTINE OLDSOLUTION

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  include 'parallel.h'

  !---------------------------------------------------------------------------------                
  ! Local variables                                                                                 
  INTEGER ::  i, j, k,n
  !---------------------------------------------------------------------------------                

!----STORE OLD SOLUTION
  do i = -1,drange(2)+1
     do k = -1,drange(4)+1
        do j = 0,drange(6)
           do n = 1,neqmax
              oldsoln(i,k,j,n) = f(i,k,j,n)
              dff(i,k,j,n) = zero
           enddo
           do n = 1,ndim
              oldsoln(i,k,j,n+maxsize) = q(i,k,j,n)
           enddo
        enddo
     enddo
  enddo

RETURN
END SUBROUTINE OLDSOLUTION
! *********************************************
