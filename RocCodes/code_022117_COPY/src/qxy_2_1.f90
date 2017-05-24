! ********************************************************************
!     SUBROUTINE UPDATE_QXY
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
! 
! ********************************************************************
SUBROUTINE UPDATE_QXY

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  include "parallel.h"

!--------------------------------------------------------------------
! Local Variables
  INTEGER :: i, j, k, m,pik(4)
  INTEGER :: eqn,maxcount, col(4)
  REAL*8 :: res_out(3)
  INTEGER :: eqn_loop(3)
!--------------------------------------------------------------------
  eqn_loop = (/4,3,4/)
  res_out = zero

  DO M=1,NDIM 

    if(eqn_loop(m) <=  0) then
       q(:,:,:,m) = 0
       cycle
    endif

    eqn = m

    CALL VELOCITY(eqn)
    CALL COEFFMATRIXQ(eqn)
    call EXP_RHSQ(eqn)
    !call FourierExtVelUpdate(2,eqn)

    it_gss=0
    res_gas=10.0
    do while (abs(res_gas) >= tol_gss .AND. it_gss <= eqn_loop(m))
       it_gss = it_gss+1
       CALL LGAUSSQ(eqn,it_gss)
    enddo
    res_out(eqn) = abs(res_gas)
  ENDDO
  maxres(4) = maxval(res_out(1:3))
  
  if(myid == 0.and. mod(ncyc_run,writemod) == 0) &
       &write(*,'(i6,a,1p10e12.4)')ncyc,'Q residuals',res_out(1:3)
  

!-------------------------------------------------------------------------  
  RETURN
END SUBROUTINE UPDATE_QXY
!*************************************************************************


!*************************************************************************
SUBROUTINE COEFFMATRIXQ(eqn)  
USE GLOBAL_DATA

  IMPLICIT NONE

!-------------------------------------------------------------------------  
! Local variables
  INTEGER ::  i, j, k, eqn, ys,ye,xs,xe,zs,ze,L
  REAL*8  :: dx1, dx2, dy1, dy2, dz1, dz2
  REAL*8  :: dx1a,dxcp1a,dxsa,dx4a,dz1a,dzcp1a,dzsa,dz4a,dy1a,dy2a,dysa,dy4a
  REAL*8  :: dxy4,dzy4,dycp1a,dyya,dyyam,ey_pointv
  REAL*8  :: cp_prndtl,dtcpp,coe
  REAL*8  :: fxbar,fzbar,fzavg,fzavgm,fxavg,fxavgm
  REAL*8  :: c_elim,cjswall,cj1wall,cj2wall,dterm,dterm1,c1,sumc_abs,abv
  LOGICAL :: drop2first(3)

  INTEGER :: ishift,kshift,othereqn,sldeqn(8),m,ee
  REAL*8,DIMENSION(2)  :: fbar,fbarinc,favginc,favgminc,Rdxdiffusion
!-------------------------------------------------------------------------  

  cp_prndtl=cp/pr

  ys = drange(5)
  ye = drange(6)
  if(dSerialRange(5) == 0) ys = ys + 1
  if(dSerialRange(6) == ny) ye = ye - 1
  zs = drange(3)
  ze = drange(4)
  xs = drange(1)
  xe = drange(2)

! Initialize variables
  dx1 = dx
  dy1 = dy
  dz1 = dz
  
!.........................................................................
! cp_prndtl is included into dx2a dz2a dy2a, used for convection terms only
!.........................................................................
  dx1a = one/dx1
  dxcp1a = cp_prndtl*half*dx1a
  dxsa = dx1a*dx1a
  dx4a = quarter*dx1a
  dy1a = one/dy1 
  dy2a = cp_prndtl*half*dy1a
  dysa = dy1a*dy1a
  dy4a = quarter*dy1a
  dz1a = one/dz1 
  dzcp1a = cp_prndtl*half*dz1a
  dzsa = dz1a*dz1a
  dz4a = quarter*dz1a

  cjswall = 8.0d0/3.0d0
  cj1wall = two
  cj2wall = third

  if (eqn <= 2) then

     mcomp  = 1;
     mgradj = 2
     mgradjM = 0
     ishift = merge(1,0,eqn == 1)
     kshift = merge(1,0,eqn == 2)
     Rdxdiffusion = (/dxsa,dzsa/)
     othereqn = merge(1,2,eqn == 2)
     if(eqn == 1) then
        sldeqn = (/ (i,i=1,8)/)
     elseif(eqn == 2) then
        sldeqn = (/ (i,i=5,8),(i,i=1,4)/)
     endif

  DO j = ys, ye
     dxy4 = dx1a*dy4a*detady(j)
     dzy4 = dz1a*dy4a*detady(j)
     dycp1a = dy2a*detady(j)
     Dyya =  dysa * detadya(j) * detady(j)
     Dyyam = dysa * detadya(j-1) * detady(j)

     do k = zs-1, ze+1
        do i = xs-1, xe+1
           MQterm(i,k,0) = (one +  MQchipr(j-1,2)*(MQphi(i+ishift,k+kshift)+MQphi(i,k))*half)
           do m = 1,2
              MQterm(i,k,m) = (one + MQchipr(j,m)*(MQphi(i+ishift,k+kshift)+MQphi(i,k))*half)
           enddo
        enddo
     enddo

  DO k = zs, ze
  DO i = xs, xe
!
     quartuple =(/myid,i,k,j/)
     qpick = (/28,3,3,1/)
!..................................................................................
! COEFF(EQN,1,...) is on the LHS, positive sign, and it is the inverse of the
! matrix coefficient of the i,k,j point
!..................................................................................
     dtcpp=dtx(i,k,j)/cp_prndtl/MQterm(i,k,mcomp)

     if(eqn == 1) then
        fbar = (/ dphidxa(i,k),half*(dphidz(i+1,k)+dphidz(i,k)) /)
        fbarinc = fbar * (/dxy4,dzy4 /)
        favginc = (/dphidx(i+1,k),0.5*(dphidza(i+1,k) + dphidza(i,k))/) * (/dxy4,dzy4 /)
        favgminc = (/dphidx(i,k),0.5*(dphidza(i+1,k-1) + dphidza(i,k-1))/) * (/dxy4,dzy4 /)
     elseif(eqn == 2) then
        fbar = (/half*(dphidx(i,k+1)+dphidx(i,k)),dphidza(i,k) /)
        fbarinc = fbar * (/dxy4,dzy4 /)
        favginc = (/0.5*(dphidxa(i,k+1) + dphidxa(i,k)),dphidz(i,k+1)/) * (/dxy4,dzy4 /)
        favgminc = (/0.5*(dphidxa(i-1,k+1) + dphidxa(i-1,k)),dphidz(i,k)/) * (/dxy4,dzy4 /)
     endif
     dyterm(1) = (one + MQchi(j-1,mgradj)*MQchi(j-1,mcomp)*(sum(fbar**2) ))&
          &/MQterm(i,k,mgradjM)*Dyyam
     dyterm(2) = (one + MQchi(j,mgradj)  *MQchi(j+1,mcomp)*(sum(fbar**2) ))&
          &/MQterm(i,k,mgradj)*Dyya

     convTerm(1,1) = (lbdavg(i,k,j,eqn+1)*MQchi(j,mgradj)*MQterm(i-1,k,mgradj)/MQterm(i,k,mgradj)&
          &-lbdavg(i,k,j-1,eqn+1)*MQchi(j-1,mgradj)*MQterm(i-1,k,mgradjM)/MQterm(i, k, mgradjM))&
          &*fbarinc(1)
     convTerm(1,2) = (lbdavg(i,k,j,eqn+1)*MQchi(j,mgradj)*MQterm(i+1,k,mgradj)/MQterm(i,k,mgradj)&
          &-lbdavg(i,k,j-1,eqn+1)*MQchi(j-1,mgradj)*MQterm(i+1,k,mgradjM)/MQterm(i, k, mgradjM))&
          &*fbarinc(1)
     convTerm(2,1) = (lbdavg(i,k,j,eqn+1)*MQchi(j,mgradj)*MQterm(i,k-1,mgradj)/MQterm(i,k,mgradj)&
          &-lbdavg(i,k,j-1,eqn+1)*MQchi(j-1,mgradj)*MQterm(i,k-1,mgradjM)/MQterm(i, k, mgradjM))&
          &*fbarinc(2)
     convTerm(2,2) = (lbdavg(i,k,j,eqn+1)*MQchi(j,mgradj)*MQterm(i,k+1,mgradj)/MQterm(i,k,mgradj)&
          &-lbdavg(i,k,j-1,eqn+1)*MQchi(j-1,mgradj)*MQterm(i,k+1,mgradjM)/MQterm(i, k, mgradjM))&
          &*fbarinc(2)
     convTerm(3,1) = (lambdag(i+ishift,k+kshift,j)*favginc(eqn)-lambdag(i,k,j)*favgminc(eqn)+ &
          lbdavg(i,k,j,1)*favginc(othereqn)-lbdavg(i-kshift,k-ishift,j,1)*favgminc(othereqn) )&
          &*MQchi(j-1,mcomp)
     convTerm(3,2) = (lambdag(i+ishift,k+kshift,j)*favginc(eqn)-lambdag(i,k,j)*favgminc(eqn)+ &
          lbdavg(i,k,j,1)*favginc(othereqn)-lbdavg(i-kshift,k-ishift,j,1)*favgminc(othereqn) )&
          &*MQchi(j+1,mcomp)


     if(skipDeformationFluid) then
        convTerm = zero
        dyterm = one
     endif

     uwvterm(1) = (uvel(i,k,j))*dxcp1a*MQterm(i,k,mcomp)  !the last term is because of dtcpp
     uwvterm(2) = (wvel(i,k,j))*dzcp1a*MQterm(i,k,mcomp)
     uwvterm(3) = (vvel(i,k,j))*dycp1a*MQterm(i,k,mcomp)

     diffterm(eqn,1) =  MQterm(i-ishift,k-kshift,mcomp)*lambdag(i,k,j)*Rdxdiffusion(eqn)
     diffterm(eqn,2) =  MQterm(i+ishift,k+kshift,mcomp)*lambdag(i+ishift,k+kshift,j)&
          &*Rdxdiffusion(eqn)
     diffterm(othereqn,1) =  MQterm(i-kshift,k-ishift,mcomp)*lbdavg(i-kshift,k-ishift,j,1)&
          &*Rdxdiffusion(othereqn)
     diffterm(othereqn,2) =  MQterm(i+kshift,k+ishift,mcomp)*lbdavg(i,k,j,1)*Rdxdiffusion(othereqn)
     diffterm(3,1) =  lbdavg(i,k,j-1,eqn+1)*dyterm(1)
     diffterm(3,2) =  lbdavg(i,k,j,eqn+1)  *dyterm(2)

     coeff(eqn,2,i,k,j) = uwvterm(1) + convTerm(1,1) + diffterm(1,1)
     coeff(eqn,3,i,k,j) =-uwvterm(1) - convTerm(1,2) + diffterm(1,2)
                     
     drop2first(1) = any(coeff(eqn,2:3,i,k,j)< zero)

     coeff(eqn,4,i,k,j) = uwvterm(2) + convTerm(2,1) + diffterm(2,1)
     coeff(eqn,5,i,k,j) =-uwvterm(2) - convTerm(2,2) + diffterm(2,2)
                      
     drop2first(2) = any(coeff(eqn,4:5,i,k,j)< zero)

     coeff(eqn,6,i,k,j) = uwvterm(3) + convTerm(3,1) + diffterm(3,1)
     coeff(eqn,7,i,k,j) =-uwvterm(3) - convTerm(3,2) + diffterm(3,2)
                     
     drop2first(3) = any(coeff(eqn,6:7,i,k,j)< zero)
     
     iconv(i,k,j,1:ndim) = merge( (/ (1,L=1,ndim)/),(/ (0,L=1,ndim)/),uwvterm(1:ndim)<=0)
     dconv(i,k,j,1:ndim) = zero
     if(drop2first(1) ) then   !drop to first order   
        dconv(i,k,j,1) = abs(uwvterm(1))
        coeff(eqn,2:3,i,k,j) =  coeff(eqn,2:3,i,k,j) + dconv(i,k,j,1)
     endif
     if(drop2first(3) ) then   !drop to first order  
        dconv(i,k,j,3) = abs(uwvterm(3))
        coeff(eqn,6:7,i,k,j) =  coeff(eqn,6:7,i,k,j) + dconv(i,k,j,3)
     endif
     if(drop2first(2) ) then   !drop to first order                                           
        dconv(i,k,j,2) = abs(uwvterm(2))
        coeff(eqn,4:5,i,k,j) =  coeff(eqn,4:5,i,k,j) + dconv(i,k,j,2)
     endif
     dconv(i,k,j,:) = dconv(i,k,j,:)*dtcpp
!!$     endif

     if((any(coeff(eqn,2:7,i,k,j) < zero) .and. j <= JRFLO(eqn)) ) then
        write(*,*) ''
        write(*,'(a,i3,3i3,a,i2,3L2,a,1p7e12.4,a,1p3e12.4,a,1p3e12.4,a,&
             &1p3e12.4,a,1p3e12.4,a,1p4e12.4,a,1p4e12.4,a,1p14e12.4)')&
             &'Myid is',myid,i,k,j,'EQN',eqn,drop2first(1:3),&
             &'Coeffs(1-7) <0',coeff(eqn,1:7,i,k,j)
     endif

     coeff(eqn,2:7,i,k,j) = coeff(eqn,2:7,i,k,j)*dtcpp


     coe = dtcpp;ee = eqn;
     csld(sldeqn(1),i,k,j) = -(MQchi(j-1,mcomp)*lambdag(i,k,j)*favgminc(ee) &
          &+lbdavg(i,k,j-1,ee+1)*fbarinc(ee)*MQchi(j-1,mgradj)*MQterm(i-ishift,k-kshift,mgradjM)&
          &/MQterm(i,k,mgradjM))*coe
     csld(sldeqn(2),i,k,j) =  (MQchi(j-1, mcomp)*lambdag(i+ishift,k+kshift,j)*favginc(ee) &
          &+lbdavg(i,k,j-1,ee+1)*fbarinc(ee)*MQchi(j-1,mgradj)*MQterm(i+ishift,k+kshift,mgradjM)&
          &/MQterm(i,k,mgradjM))*coe
     csld(sldeqn(3),i,k,j) = -(MQchi(j+1,mcomp)*lambdag(i+ishift,k+kshift,j)*favginc(ee) &
          &+lbdavg(i,k,j,ee+1)*fbarinc(ee)*MQchi(j,mgradj)*MQterm(i+ishift,k+kshift,mgradj)&
          &/MQterm(i,k,mgradj))*coe
     csld(sldeqn(4),i,k,j) =  (MQchi(j+1, mcomp)*lambdag(i, k, j)*favgminc(ee) &
          &+lbdavg(i,k,j,ee+1)*fbarinc(ee)*MQchi(j,mgradj)*MQterm(i-ishift,k-kshift,mgradj)&
          &/MQterm(i,k,mgradj))*coe

     coe = dtcpp;ee = othereqn;
     csld(sldeqn(5),i,k,j) = -(MQchi(j-1,mcomp)*lbdavg(i-kshift,k-ishift,j,1)*favgminc(ee) &
          &+lbdavg(i,k,j-1,eqn+1)*fbarinc(ee)*MQchi(j-1,mgradj)*MQterm(i-kshift,k-ishift,mgradjM)&
          &/MQterm(i,k,mgradjM))*coe
     csld(sldeqn(6),i,k,j) =  (MQchi(j-1,mcomp)*lbdavg(i,k,j,1)*favginc(ee) &
          &+lbdavg(i,k,j-1,eqn+1)*fbarinc(ee)*MQchi(j-1,mgradj)*MQterm(i+kshift,k+ishift,mgradjM)&
          &/MQterm(i,k,mgradjM))*coe
     csld(sldeqn(7),i,k,j) = -(MQchi(j+1,mcomp)*lbdavg(i,k,j,1)*favginc(ee) &
          &+ lbdavg(i,k,j,eqn+1) *fbarinc(ee)*MQchi(j,mgradj)*MQterm(i+kshift,k+ishift,mgradj)&
          &/MQterm(i,k,mgradj))*coe
     csld(sldeqn(8),i,k,j) = (MQchi(j+1,mcomp)*lbdavg(i-kshift,k-ishift,j,1)*favgminc(ee)&
          &+ lbdavg(i,k,j,eqn+1) *fbarinc(ee)*MQchi(j,mgradj)*MQterm(i-kshift,k-ishift,mgradj)&
          &/MQterm(i,k,mgradj))*coe

     coeff(eqn,1,i,k,j) = one +  sum(coeff(eqn,2:7,i,k,j))   + sum(csld(1:8,i,k,j))

     if(j == 1) then
        c_elim = coeff(eqn,6,i,k,j)
        coeff(eqn,7,i,k,j) = coeff(eqn,7,i,k,j) + cj2wall*c_elim
        coeff(eqn,1,i,k,j) = coeff(eqn,1,i,k,j) + cj1wall*c_elim     !cj1wall is the negative of the coefficient of u1 in u0
        dqdt(i,k,j,eqn) = dqdt(i,k,j,eqn) + cjswall*c_elim*srfqbnd(i,k,eqn)
        coeff(eqn,6,i,k,j) = zero

        c_elim = csld(1,i,k,j)
        csld(4,i,k,j) = csld(4,i,k,j) + cj2wall*c_elim
        coeff(eqn,2,i,k,j) = coeff(eqn,2,i,k,j) - cj1wall*c_elim
        dqdt(i,k,j,eqn) = dqdt(i,k,j,eqn) + cjswall*c_elim*srfqbnd(i-1,k,eqn)
        csld(1,i,k,j) = zero
        
        c_elim = csld(2,i,k,j)
        csld(3,i,k,j) = csld(3,i,k,j) + cj2wall*c_elim
        coeff(eqn,3,i,k,j) = coeff(eqn,3,i,k,j) - cj1wall*c_elim
        dqdt(i,k,j,eqn) = dqdt(i,k,j,eqn) + cjswall*c_elim*srfqbnd(i+1,k,eqn)
        csld(2,i,k,j) = zero
        
        c_elim = csld(5,i,k,j)
        csld(8,i,k,j) = csld(8,i,k,j) + cj2wall*c_elim
        coeff(eqn,4,i,k,j) = coeff(eqn,4,i,k,j) - cj1wall*c_elim
        dqdt(i,k,j,eqn) = dqdt(i,k,j,eqn) + cjswall*c_elim*srfqbnd(i,k-1,eqn)
        csld(5,i,k,j) = zero
        
        c_elim = csld(6,i,k,j)
        csld(7,i,k,j) = csld(7,i,k,j) + cj2wall*c_elim
        coeff(eqn,5,i,k,j) = coeff(eqn,5,i,k,j) - cj1wall*c_elim
        dqdt(i,k,j,eqn) = dqdt(i,k,j,eqn) + cjswall*c_elim*srfqbnd(i,k+1,eqn)
        csld(6,i,k,j) = zero
     endif


  ENDDO
  ENDDO
  ENDDO   
!
  elseif (eqn == 3) then
     mcomp = 2;
     mgradj = 1
     mgradjM = 0

  DO j = ys, ye
     ey_pointv = detadya(j)
     dycp1a = dy2a*ey_pointv
     Dyya =  dysa*detady(j+1)*ey_pointv
     Dyyam = dysa*detady(j)*ey_pointv
     dxy4 = dx1a*dy4a*ey_pointv
     dzy4 = dz1a*dy4a*ey_pointv
     do k = zs-1, ze+1
        do i = xs-1, xe+1
           MQterm(i,k,mcomp)   = (one +  MQchipr(j,2)*  MQphi(i,k))
           MQterm(i,k,mgradj)  = (one +  MQchipr(j+1,1)*MQphi(i,k))
           MQterm(i,k,mgradjM) = (one +  MQchipr(j,1)*  MQphi(i,k))
        enddo
     enddo

  DO k = zs, ze
  DO i = xs, xe

     dtcpp=dtx(i,k,j)/cp_prndtl/MQterm(i,k,mcomp)

     fbar = (/ dphidx(i,k), dphidz(i,k)/)
     fbarinc = fbar * (/dxy4,dzy4 /)
     dyterm(1) = (one + MQchi(j,mgradj)  *MQchi(j-1,mcomp)*(sum(fbar**2) ))&
          &/MQterm(i,k,mgradjM)*Dyyam
     dyterm(2) = (one + MQchi(j+1,mgradj)*MQchi(j+1,mcomp)*(sum(fbar**2) ))&
          &/MQterm(i,k,mgradj )*Dyya


     convTerm(1,1) = (lambdag(i,k,j+1)*MQchi(j+1,mgradj)*MQterm(i-1,k,mgradj)/MQterm(i,k,mgradj)&
          &-lambdag(i,k,j)*MQchi(j, mgradj)*MQterm(i-1,k,mgradjM)/MQterm(i,k,mgradjM))&
          &*fbarinc(1)
     convTerm(1,2) = (lambdag(i,k,j+1)*MQchi(j+1,mgradj)*MQterm(i+1,k,mgradj)/MQterm(i,k,mgradj)&
          &-lambdag(i,k,j)*MQchi(j,mgradj)*MQterm(i+1,k,mgradjM)/MQterm(i,k,mgradjM))&
          &*fbarinc(1)
     convTerm(2,1) = (lambdag(i,k,j+1)*MQchi(j+1,mgradj)*MQterm(i,k-1,mgradj)/MQterm(i,k,mgradj)&
          &-lambdag(i,k,j)*MQchi(j, mgradj)*MQterm(i,k-1,mgradjM)/MQterm(i,k,mgradjM))&
          &*fbarinc(2)
     convTerm(2,2) = (lambdag(i,k,j+1)*MQchi(j+1,mgradj)*MQterm(i,k+1,mgradj)/MQterm(i,k,mgradj)&
          &-lambdag(i,k,j)*MQchi(j,mgradj)*MQterm(i,k+1,mgradjM)/MQterm(i,k,mgradjM))&
          &*fbarinc(1)
     convTerm(3,1) = MQchi(j-1, mcomp)*( &
          &(lbdavg(i,k,j,3)*dphidza(i,k)-lbdavg(i,k-1,j,3)*dphidza(i,k-1))*dzy4  &
          +(lbdavg(i,k,j,2)*dphidxa(i,k)-lbdavg(i-1,k,j,2)*dphidxa(i-1,k))*dxy4)
     convTerm(3,2) = MQchi(j+1, mcomp)*( &
          &(lbdavg(i,k,j,3)*dphidza(i,k)-lbdavg(i,k-1,j,3)*dphidza(i,k-1))*dzy4  &
          +(lbdavg(i,k,j,2)*dphidxa(i,k)-lbdavg(i-1,k,j,2)*dphidxa(i-1,k))*dxy4)

     uwvterm(1) = (uvel(i,k,j))*dxcp1a
     uwvterm(2) = (wvel(i,k,j))*dzcp1a
     uwvterm(3) = (vvel(i,k,j))*dycp1a


     diffterm(1,1) =  lbdavg(i-1,k,j,2)*dxsa*MQterm(i-1,k,mcomp)
     diffterm(1,2) =  lbdavg(i,  k,j,2)*dxsa*MQterm(i+1,k,mcomp)
     diffterm(2,1) =  lbdavg(i,k-1,j,3)*dzsa*MQterm(i,k-1,mcomp)
     diffterm(2,2) =  lbdavg(i,  k,j,3)*dzsa*MQterm(i,k+1,mcomp)
     diffterm(3,1) =  lambdag(i,k,j  )*dyterm(1)
     diffterm(3,2) =  lambdag(i,k,j+1)*dyterm(2)

     if(skipDeformationFluid) then
        convTerm = zero
        dyterm = one
     endif

     coeff(eqn,2,i,k,j) = uwvterm(1) + convTerm(1,1) + diffterm(1,1)
     coeff(eqn,3,i,k,j) =-uwvterm(1) - convTerm(1,2) + diffterm(1,2)
                     
     drop2first(1) = any(coeff(eqn,2:3,i,k,j)< zero)

     coeff(eqn,4,i,k,j) = uwvterm(2) + convTerm(2,1) + diffterm(2,1)
     coeff(eqn,5,i,k,j) =-uwvterm(2) - convTerm(2,2) + diffterm(2,2)
                      
     drop2first(2) = any(coeff(eqn,4:5,i,k,j)< zero)

     coeff(eqn,6,i,k,j) = uwvterm(3) + convTerm(3,1) + diffterm(3,1)
     coeff(eqn,7,i,k,j) =-uwvterm(3) - convTerm(3,2) + diffterm(3,2)
                     
     drop2first(3) = any(coeff(eqn,6:7,i,k,j)< zero)

     iconv(i,k,j,1:ndim) = merge( (/ (1,L=1,ndim)/),(/ (0,L=1,ndim)/),uwvterm(1:ndim)<=0)
     dconv(i,k,j,1:ndim) = zero
     if(drop2first(1) ) then   !drop to first order 
        dconv(i,k,j,1) = abs(uwvterm(1))
        coeff(eqn,2:3,i,k,j) =  coeff(eqn,2:3,i,k,j) + dconv(i,k,j,1)
     endif
     if(drop2first(3) ) then   !drop to first order   
        dconv(i,k,j,3) = abs(uwvterm(3))
        coeff(eqn,6:7,i,k,j) =  coeff(eqn,6:7,i,k,j) + dconv(i,k,j,3)
     endif
     if(drop2first(2) ) then   !drop to first order                                           
        dconv(i,k,j,2) = abs(uwvterm(2))
        coeff(eqn,4:5,i,k,j) =  coeff(eqn,4:5,i,k,j) + dconv(i,k,j,2)
     endif
     dconv(i,k,j,:) = dconv(i,k,j,:)*dtcpp


     if (j == JRFLO(3)-1)  then
        coeff(eqn,6,i,k,j)  = sum(convTerm(3,1:2) + diffterm(3,1:2))
        coeff(eqn,7,i,k,j)  = zero
     endif
     if (j == ny-1) dconv(i,k,j,3) = zero


     if((any(coeff(eqn,2:7,i,k,j) < zero) .and. j <= JRFLO(eqn)) ) then
        write(*,*) ''
        write(*,'(a,i3,3i3,a,i2,3L2,a,1p7e12.4,a,1p3e12.4,a,1p3e12.4,a,&
             &1p3e12.4,a,1p3e12.4,a,1p4e12.4,a,1p4e12.4,a,1p14e12.4)')&
             &'Myid is',myid,i,k,j,'EQN',eqn,drop2first(1:3),&
             &'Coeffs(1-7) <0',coeff(eqn,1:7,i,k,j)
     endif

     coeff(eqn,2:7,i,k,j) = coeff(eqn,2:7,i,k,j)*dtcpp


!..................................................................................     
! TO save memory use csld to store the orhter terms (corners) in the jacobian Matrix
!..................................................................................     
     coe = dxy4*dtcpp
     csld(1,i,k,j) = -(lambdag(i,k,j)*dphidx(i,k)*&
          &MQchi(j,mgradj)*MQterm(i-1,k,mgradjM)/MQterm(i,k,mgradjM)&
          &+lbdavg(i-1,k,j,2)*dphidxa(i-1,k)*MQchi(j-1,mcomp))*coe

     csld(2,i,k,j) =  (lambdag(i,k,j)*dphidx(i, k)*&
          &MQchi(j,mgradj)*MQterm(i+1, k,mgradjM)/MQterm(i,k,mgradjM)&
          &+lbdavg(i,k,j,2)*dphidxa(i,k)*MQchi(j-1,mcomp)  )*coe

     csld(3,i,k,j) = -(lambdag(i,k,j+1)*dphidx(i, k)*&
          &MQchi(j + 1, mgradj)*MQterm(i + 1, k, mgradj)/MQterm(i, k, mgradj)&
          &+lbdavg(i,k,j,2)*dphidxa(i,k)*MQchi(j+1,mcomp))*coe
     csld(4,i,k,j) = (lambdag(i,k,j+1)*dphidx(i, k)*&
          &MQchi(j + 1, mgradj)*MQterm(i - 1, k, mgradj)/MQterm(i, k, mgradj)&
          &+lbdavg(i-1,k,j,2)*dphidxa(i-1,k)*MQchi(j+1,mcomp))*coe

     coe = dzy4*dtcpp
     csld(5,i,k,j) = -(lambdag(i,k,j)*dphidz(i,k)*&
          &MQchi(j,mgradj)*MQterm(i,k-1,mgradjM)/MQterm(i,k,mgradjM)&
          &+lbdavg(i,k-1,j,3)*dphidza(i,k-1)*MQchi(j-1,mcomp))*coe

     csld(6,i,k,j) =  (lambdag(i,k,j)*dphidz(i,k)*&
          &MQchi(j,mgradj)*MQterm(i,k+1,mgradjM)/MQterm(i,k,mgradjM)&
          &+lbdavg(i,k,j,3)*dphidza(i,k)*MQchi(j-1,mcomp) )*coe

     csld(7,i,k,j) = -(lambdag(i,k,j+1)*dphidz(i,k)*&
          &MQchi(j+1,mgradj)*MQterm(i,k+1,mgradj)/MQterm(i,k,mgradj)&
          &+lbdavg(i,k,j,3)*dphidza(i,k)*MQchi(j+1,mcomp))*coe
     csld(8,i,k,j) = (lambdag(i,k,j+1)*dphidz(i,k)*&
          &MQchi(j+1,mgradj)*MQterm(i,k-1,mgradj)/MQterm(i,k,mgradj)&
          &+lbdavg(i,k-1,j,3)*dphidza(i,k-1)*MQchi(j+1,mcomp))*coe

     coeff(eqn,1,i,k,j) = one +  sum(coeff(eqn,2:7,i,k,j))   + sum(csld(1:8,i,k,j))

  ENDDO
  ENDDO
  ENDDO   

  endif  !if (eqn == 1)

  if(skipDeformationFluid) csld = zero
 
!--------------------------------------------------------------------
   RETURN
END SUBROUTINE COEFFMATRIXQ
!*********************************************************************


!*********************************************************************
SUBROUTINE LGAUSSQ(eqn,icount)

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

  include "parallel.h"

  !---------------------------------------------------------------------
  ! Dummy variables
  INTEGER, INTENT(IN) :: eqn
  ! Local variables
  ! c_i::center point jacobian, b_i and d_i::corners
  INTEGER ::  i, j, k, jp,jp1, js, jmax,idir,icount

  INTEGER ::   ierrwg,m ,wunit
  INTEGER ::  xs,xe,zs,ze,ys,ye,indx,icolor, col(7)
  INTEGER,DIMENSION(3) :: ikj,ikj_P, ikj_P1,ikj_M1,ikj_M2,ikjP,ikjM
  REAL*8  ::  coe,c2,c3,c4,c5,b1,b2,b3,b4,d1,d2,d3,d4
  REAL*8  ::  oldq,rhsg,resid,omg_uvw,vv(12),cc(12)
  REAL*8  ::  s1,s2,s3,sig1,sig2,corrflux,nuy
  logical :: outcome
  !---------------------------------------------------------------------

  wunit = 6
  ! Initialize variables
  ys = drange(5)
  ye = drange(6)
  if(dSerialRange(5).eq.0) ys = ys + 1
  if(dSerialRange(6).eq.ny) ye = ye - 1
  col(3:4)=eqn
  col(5:6) = drange(5:6)
  col(7) = 0
  omg_uvw = 1.00d0

  resid= -100.d0

  if(JRFLO(eqn) >= ny ) then
     jmax = ye+1
     divt(:,:,jmax) = 0d0
  else
     jmax = JRFLO(eqn)
     if(eqn == 1) then
        q(:,:,jmax,eqn) = Uimposed(2)  ! new time step value
     elseif(eqn == 2) then
        q(:,:,jmax,eqn) = zero
     else
        jmax = JRFLO(eqn)
        divt(:,:,jmax) = 0d0   !Neuman boundary conditions
     endif
  endif
  js = ys-1

  do icolor = 1,2
     col(1:2) = icolor
     do indx = 1 , finest_mesh%nclr(icolor)
        i =  finest_mesh%iclr(indx,icolor)
        k =  finest_mesh%kclr(indx,icolor)

        !--- setup Thomas algorithm
        fy(js) = newqbnd(i,k,eqn)
        ay(js) = 0.0d0
        cy(js) = 0.0d0
        by(js) = 1.0d0
        do j = ys, jmax-1
           jp1 = merge(j+1,j-1, j<jmax-1 .OR. eqn<3)

           ay(j) = - coeff(eqn,6,i,k,j)
           cy(j) = - coeff(eqn,7,i,k,j)
           by(j) = coeff(eqn,1,i,k,j)
           !
           ikj = (/i,k,j/)
           corrflux = zero
           if(.not. skipTVD) then
              do idir = 1,ndim
                 nuy = dconv(i,k,j,idir)

                 if(nuy < 1.d-9) CYCLE

                 ikj_P = ikj;ikj_P(idir) = ikj(idir) + iconv(i,k,j,idir)
                 ikj_P1 = ikj_P; ikj_P1(idir) = ikj_P(idir)+1
                 ikj_M1 = ikj_P; ikj_M1(idir) = ikj_P(idir)-1
                 ikj_M2 = ikj_P; ikj_M2(idir) = ikj_P(idir)-2

                 s1 = (q(ikj_P1(1),ikj_P1(2),ikj_P1(3),eqn) - q(ikj_P(1),ikj_P(2),ikj_P(3),eqn))
                 s2 = (q(ikj_P(1),ikj_P(2),ikj_P(3),eqn) - q(ikj_M1(1),ikj_M1(2),ikj_M1(3),eqn))
                 s3 = (q(ikj_M1(1),ikj_M1(2),ikj_M1(3),eqn) - q(ikj_M2(1),ikj_M2(2),ikj_M2(3),eqn))
                 sig1 = (sign(half,s1)+sign(half,s2))*min(abs(s1),abs(s2))
                 sig2 = (sign(half,s2)+sign(half,s3))*min(abs(s2),abs(s3))

                 corrflux =  corrflux + nuy * (sig1-sig2) 
                 if(ncyc_run < ncyc_debug) then
                    rate(i,k,j,idir) = nuy * (sig1-sig2)
                 endif
              enddo
           endif

           vv = (/q(i-1,k,j,eqn),q(i+1,k,j,eqn),q(i,k-1,j,eqn),q(i,k+1,j,eqn),&
                &q(i-1,k,j-1,eqn),q(i+1,k,j-1,eqn),q(i+1,k,jp1,eqn),q(i-1,k,jp1,eqn),&
                &q(i,k-1,j-1,eqn),q(i,k+1,j-1,eqn),q(i,k+1,jp1,eqn),q(i,k-1,jp1,eqn) /)
           cc(1:4) = coeff(eqn,2:5,i,k,j);cc(5:12) = csld(1:8,i,k,j)
!
           fy(j)= sum(cc*vv) + dqdt(i,k,j,eqn) -  corrflux

           oldq = fy(j) - ay(j)*q(i,k,j-1,eqn) - cy(j)*q(i,k,j+1,eqn) - by(j)*q(i,k,j,eqn)

           quartuple =(/myid,i,k,j/)
           qpick = (/4, 2, 10, 1/)

           resid = max(resid,abs(oldq))
        end do
        if(eqn < 3) then
           fy(jmax) = q(i,k,jmax,eqn)
           ay(jmax) = 0.0d0
           cy(jmax) = 0.0d0
           by(jmax) = 1.0d0
        else
           fy(jmax) = divt(i,k,jmax)
           ay(jmax) = -1.0d0
           cy(jmax) = 0.0d0
           by(jmax) = 1.0d0
        endif
        
!----   THOMAS algorithm
        do j=js+1,jmax
           ay(j) = ay(j)/by(j-1)
           by(j) = by(j)-ay(j)*cy(j-1)
           fy(j) = fy(j)-ay(j)*fy(j-1)
        enddo

        fy(jmax) = fy(jmax)/by(jmax)
        q(i,k,jmax,eqn) = fy(jmax)

        do j=jmax-1,js,-1
           fy(j) = (fy(j)-cy(j)*fy(j+1))/by(j)
           q(i,k,j,eqn) = fy(j)          

           if(ncyc_run < ncyc_debug) then
              call NAN_VAL(q(i,k,j,eqn),outcome,1)
              if(outcome) then
                 write(wunit,*)'============================================'
                 write(wunit,*)'DEBUG SECTION FOR VELOCITY SOLVER,eqn =',eqn
                 write(wunit,*)'N_CYCLE',ncyc
                 write(wunit,*)'Myid',myid
                 write(wunit,*)'ITER',icount,icolor
                 write(wunit,'(a,1p3e12.4,a,2l2)')'Time order:coe_exp,time_coeff,dtx(i,k,j)',&
                      &coe_exp,time_coeff,dtx(i,k,j),'Flags{skiptvd,is_first_order}',skipTVD,is_firstorder
                 write(wunit,*)'POS REL This proc',i,k,j
                 write(wunit,*)'grid limits',drange(2),drange(4),drange(6)
                 write(wunit,*)'POS ABS',i+ib_mg(1),k+kb_mg(1),j
                 write(wunit,'(a,1p50e12.4)')'NEW SOLN',q(i,k,j,eqn)
                 write(wunit,'(a,1p50e12.4)')'OLD SOLN',oldsoln(i,k,j,eqn+maxsize)
                 write(wunit,'(a,1p50e12.4)')'NEIGH SOL,i-',q(i-2:i,k,j,eqn)
                 write(wunit,'(a,1p50e12.4)')'NEIGH SOL,i+',q(i:i+2,k,j,eqn)
                 write(wunit,'(a,1p50e12.4)')'NEIGH SOL,k-',q(i,k-2:k,j,eqn)
                 write(wunit,'(a,1p50e12.4)')'NEIGH SOL,k+',q(i,k:k+2,j,eqn)
                 write(wunit,'(a,1p50e12.4)')'Explicit', dqdt(i,k,j,eqn)
                 write(wunit,'(a,1p50e12.4)')'Corrective Flux', rate(i,k,j,1:ndim)
                 write(wunit,'(a,1p50e12.4)')'Corrective Velox', dconv(i,k,j,1:ndim)
                 write(wunit,'(a,1p50e12.4)')'Corner coeffs',csld(1:8,i,k,j)
                 write(wunit,'(a,1p50e12.4)')'COEFFs1-7',coeff(eqn,1:7,i,k,j)
                 write(wunit,*)'############################################'
                 call MPI_Abort(MPI_COMM_WORLD,ierrWG)
                 call MPI_finalize(MPI_COMM_WORLD)
                 stop
              endif
           endif

        enddo

     end do  !indx
     CALL parallel_swap4(col,finest_mesh)
  enddo     !icolor


  if (eqn == 3) then
     if(JRFLO(eqn) >=ny) then
        do i= drange(1)-1, drange(2)+1
           do k= drange(3)-1, drange(4)+1
              q(i,k,ny,3) = q(i,k,ny-1,3)
           enddo
        enddo
     else
        do i= drange(1)-1, drange(2)+1
           do k= drange(3)-1, drange(4)+1
              do j = JRFLO(3),ny
                 q(i,k,J,3) = q(i,k,JRFLO(3)-1,3)
              enddo
           enddo
        enddo
     endif
  endif
  
  call MPI_ALLREDUCE(resid,res_gas,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)

!---------------------------------------------------------------------
  RETURN
END SUBROUTINE LGAUSSQ
!*********************************************************************
