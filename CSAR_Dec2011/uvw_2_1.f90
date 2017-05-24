! ********************************************************************
SUBROUTINE UPDATE_UVW

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  include "parallel.h"

!--------------------------------------------------------------------
! Local Variables
  INTEGER :: i, j, k, m, n
  INTEGER :: eqn,maxcount, col(4)
  REAL*8 :: res_visc_int
!--------------------------------------------------------------------

  maxcount=merge(6,6,ncyc < ncyc_init)

  res_visc_int = zero
  DO M=1,NDIM,2               !SKIP THE SECOND EQUATION
    EQN=m
    CALL VELOCITY(eqn)
    CALL COEFFMATRIXQ(eqn)

    it_gss=0;res_gas=10.0
    do while (abs(res_gas) >= 1.d-16 .AND. it_gss <= maxcount)
       it_gss = it_gss+1
       CALL RHSQ(eqn)
    enddo

    res_visc_int = max(res_visc_int,abs(res_gas))

    CALL  DIRICHLET_BC(M,finest_mesh)

  ENDDO

 
!-------------------------------------------------------------------------  
  RETURN
END SUBROUTINE UPDATE_UVW
!*************************************************************************


!*************************************************************************
SUBROUTINE COEFFMATRIXQ(eqn)  

  USE GLOBAL_DATA

  IMPLICIT NONE

!-------------------------------------------------------------------------  
! Local variables
  INTEGER ::  i, j, eqn, k, ys,ye,xs,xe,zs,ze
  REAL*8  :: dx1, dx2, dy1, dy2, dz1, dz2
  REAL*8  :: dx1a,dx2a,dxsa,dx4a,dz1a,dz2a,dzsa,dz4a,dy1a,dy2a,dysa,dy4a
  REAL*8  :: dxy4,dzy4,dyy2,dyya,dyyam,ey_pointv
  REAL*8  :: cp_prndtl,dtcpp,coe
  REAL*8  :: fxbar,fzbar,fzavg,fzavgm,fxavg,fxavgm
  REAL*8  :: term,convx,convz,convy, uterm,vterm,wterm
  REAL*8  :: c_elim,cjswall,cj1wall,cj2wall,dterm,dterm1,c1,sumc_abs
!-------------------------------------------------------------------------  

  cp_prndtl=cp/pr

  zs = drange(3)
  ze = drange(4)
  xs = drange(1)
  xe = nxv(eqn)
  ys = drange(5)+1
  ye = nyv(eqn)

! Initialize variables
  dx1 = dx
  dy1 = dy
  dz1 = dz
  
!.........................................................................
! cp_prndtl is included into dx2a dz2a dy2a, used for convection terms only
!.........................................................................
  dx1a = one/dx1
  dx2a = cp_prndtl*half*dx1a
  dxsa = dx1a*dx1a
  dx4a = quarter*dx1a
  dy1a = one/dy1 
  dy2a = cp_prndtl*half*dy1a
  dysa = dy1a*dy1a
  dy4a = quarter*dy1a
  dz1a = one/dz1 
  dz2a = cp_prndtl*half*dz1a
  dzsa = dz1a*dz1a
  dz4a = quarter*dz1a

  cjswall=8.0d0/3.0d0;cj1wall=two;cj2wall=third;
  
  if (eqn.eq.1) then

  DO j = ys, ye
     dxy4 = dx1a*dy4a*detady(j)
     dzy4 = dz1a*dy4a*detady(j)
     dyy2 = dy2a*detady(j)
     dyya =  dysa*half*(detady(j)+detady(j+1))*detady(j)
     dyyam = dysa*half*(detady(j)+detady(j-1))*detady(j)
  DO k = zs, ze
  DO i = xs, xe
            
     dqdt(i,k,j,eqn) = dfdt(i,k,j,eqn) + rate(i,k,j,eqn)

!..................................................................................
! COEFF(1,...) is on the LHS, positive sign, and it is the inverse of the
! matrix coefficient of the i,k,j point
!..................................................................................
     dtcpp=dtx(i,k,j)/cp_prndtl
     fxbar = dphidxa(i,k)
     fzbar = 0.5*(dphidz(i+1,k)+dphidz(i,k))
     fzavg = 0.5*(dphidza(i+1,k) + dphidza(i,k))
     fzavgm= 0.5*(dphidza(i+1,k-1)+dphidza(i,k-1))
     term = one+fxbar**2+fzbar**2

     convx = (lbdavg(i,k,j,2)-lbdavg(i,k,j-1,2))*fxbar*dxy4
     convz = (lbdavg(i,k,j,2)-lbdavg(i,k,j-1,2))*fzbar*dzy4
     convy = (lambdag(i+1,k,j)*dphidx(i+1,k)-lambdag(i,k,j)*dphidx(i,k))*dxy4+ &
             (lbdavg(i,k,j,1)*fzavg-lbdavg(i,k-1,j,1)*fzavgm)*dzy4

     vterm = (vvel(i,k,j))*dyy2
     uterm = (uvel(i,k,j))*dx2a
     wterm = (wvel(i,k,j))*dz2a

     coeff(2,i,k,j) = uterm + convx + lambdag(i,k,j)*dxsa
     coeff(3,i,k,j) = -convx + lambdag(i+1,k,j)*dxsa&
                    - uterm

     coeff(4,i,k,j) = vterm + convy + lbdavg(i,k,j-1,2)*term*dyyam
     coeff(5,i,k,j) = -convy + lbdavg(i,k,j,2)*term*dyya &
                    - vterm

     coeff(6,i,k,j) = wterm + convz + lbdavg(i,k-1,j,1)*dzsa
     coeff(7,i,k,j) = -convz + lbdavg(i,k,j,1)*dzsa&
                    - wterm

     if(f1st_ord_conv) then   !drop to first order the convective terms
       coeff(2,i,k,j) =  coeff(2,i,k,j) + abs(uvel(i,k,j))*dx2a
       coeff(3,i,k,j) =  coeff(3,i,k,j) + abs(uvel(i,k,j))*dx2a
       coeff(4,i,k,j) =  coeff(4,i,k,j) + abs(vvel(i,k,j))*dyy2
       coeff(5,i,k,j) =  coeff(5,i,k,j) + abs(vvel(i,k,j))*dyy2
       coeff(6,i,k,j) =  coeff(6,i,k,j) + abs(wvel(i,k,j))*dz2a
       coeff(7,i,k,j) =  coeff(7,i,k,j) + abs(wvel(i,k,j))*dz2a
     endif

     coeff(2:7,i,k,j) = coeff(2:7,i,k,j)*dtcpp



!..................................................................................     
! TO save memory use Bm to store the orhter terms (corners) of the jacobian Matrix
!..................................................................................     
     coe = dxy4*dtcpp
     Bm(1,1,i,k,j) = -(lambdag(i,k,j)*dphidx(i,k) + lbdavg(i,k,j-1,2)*fxbar)*coe
     Bm(1,2,i,k,j) =  (lambdag(i+1,k,j)*dphidx(i+1,k) + lbdavg(i,k,j-1,2)*fxbar)*coe
     Bm(1,3,i,k,j) = -(lambdag(i+1,k,j)*dphidx(i+1,k) + lbdavg(i,k,j,2)*fxbar)*coe
     Bm(1,4,i,k,j) = - sum(Bm(1,1:3,i,k,j))

     coe = dzy4*dtcpp
     Bm(2,1,i,k,j) = -(lbdavg(i,k-1,j,1)*fzavgm + lbdavg(i,k,j-1,2)*fzbar)*coe
     Bm(2,2,i,k,j) =  (lbdavg(i,k,j,1)*fzavg + lbdavg(i,k,j-1,2)*fzbar)*coe
     Bm(2,3,i,k,j) = -(lbdavg(i,k,j,1)*fzavg + lbdavg(i,k,j,2)*fzbar)*coe
     Bm(2,4,i,k,j) = - sum(Bm(2,1:3,i,k,j))

  ENDDO
  ENDDO
  ENDDO 

  elseif (eqn.eq.2) then

  DO j = ys, ye
     dxy4 = dx1a*dy4a*detady(j)
     dzy4 = dz1a*dy4a*detady(j)
     dyy2 = dy2a*detady(j)
     dyya =  dysa*half*(detady(j)+detady(j+1))*detady(j)
     dyyam = dysa*half*(detady(j)+detady(j-1))*detady(j)
  DO k = zs, ze
  DO i = xs, xe
              
     dqdt(i,k,j,eqn) = dfdt(i,k,j,eqn) + rate(i,k,j,eqn)
   
!..................................................................................
! COEFF(1,...) is on the LHS, positive sign, and it is the inverse of the
! matrix coefficient of the i,k,j point
! cp_prndtl = cp/pr is included in dx2a as well so that when it is mult. by dtcpp is canc.
!..................................................................................
     dtcpp=dtx(i,k,j)/cp_prndtl
     fxbar = half*(dphidx(i,k+1)+dphidx(i,k))
     fzbar = dphidza(i,k)
     fxavg = half*(dphidxa(i,k+1)+dphidxa(i,k))
     fxavgm= half*(dphidxa(i-1,k+1)+dphidxa(i-1,k))
     term = one+fxbar**2+fzbar**2

     convx = (lbdavg(i,k,j,3)-lbdavg(i,k,j-1,3))*fxbar*dxy4
     convz = (lbdavg(i,k,j,3)-lbdavg(i,k,j-1,3))*fzbar*dzy4
     convy = (lambdag(i,k+1,j)*dphidz(i,k+1)-lambdag(i,k,j)*dphidz(i,k))*dzy4+ &
             (lbdavg(i,k,j,1)*fxavg-lbdavg(i-1,k,j,1)*fxavgm)*dxy4

     uterm = (uvel(i,k,j))*dx2a
     wterm = (wvel(i,k,j))*dz2a
     vterm = (vvel(i,k,j))*dyy2

     coeff(2,i,k,j) = uterm + convx + lbdavg(i-1,k,j,1)*dxsa
     coeff(3,i,k,j) = -convx + lbdavg(i,k,j,1)*dxsa&
                    - uterm

     coeff(4,i,k,j) = vterm + convy + lbdavg(i,k,j-1,3)*term*dyyam
     coeff(5,i,k,j) = - convy + lbdavg(i,k,j,3)*term*dyya&
                    - vterm

     coeff(6,i,k,j) = wterm + convz + lambdag(i,k,j)*dzsa
     coeff(7,i,k,j) = -convz + lambdag(i,k+1,j)*dzsa&
                    - wterm

     if(f1st_ord_conv) then   !drop to first order the convective terms
       coeff(2,i,k,j) =  coeff(2,i,k,j) + abs(uvel(i,k,j))*dx2a
       coeff(3,i,k,j) =  coeff(3,i,k,j) + abs(uvel(i,k,j))*dx2a
       coeff(4,i,k,j) =  coeff(4,i,k,j) + abs(vvel(i,k,j))*dyy2
       coeff(5,i,k,j) =  coeff(5,i,k,j) + abs(vvel(i,k,j))*dyy2
       coeff(6,i,k,j) =  coeff(6,i,k,j) + abs(wvel(i,k,j))*dz2a
       coeff(7,i,k,j) =  coeff(7,i,k,j) + abs(wvel(i,k,j))*dz2a
     endif

     coeff(2:7,i,k,j) = coeff(2:7,i,k,j)*dtcpp


!..................................................................................     
! TO save memory use Bm to store the orhter terms (corners) in the jacobian Matrix
!..................................................................................     
     coe = dzy4*dtcpp
     Bm(2,1,i,k,j) = -(lambdag(i,k,j)*dphidz(i,k) + lbdavg(i,k,j-1,3)*fzbar)*coe
     Bm(2,2,i,k,j) =  (lambdag(i,k+1,j)*dphidz(i,k+1) + lbdavg(i,k,j-1,3)*fzbar)*coe
     Bm(2,3,i,k,j) = -(lambdag(i,k+1,j)*dphidz(i,k+1) + lbdavg(i,k,j,3)*fzbar)*coe
     Bm(2,4,i,k,j) = - sum(Bm(2,1:3,i,k,j))

     coe = dxy4*dtcpp
     Bm(1,1,i,k,j) = -(lbdavg(i-1,k,j,1)*fxavgm + lbdavg(i,k,j-1,3)*fxbar)*coe
     Bm(1,2,i,k,j) =  (lbdavg(i,k,j,1)*fxavg + lbdavg(i,k,j-1,3)*fxbar)*coe
     Bm(1,3,i,k,j) = -(lbdavg(i,k,j,1)*fxavg + lbdavg(i,k,j,3)*fxbar)*coe
     Bm(1,4,i,k,j) = - sum(Bm(1,1:3,i,k,j))

 


  ENDDO
  ENDDO
  ENDDO   

!
  elseif (eqn.eq.3) then

  DO j = ys, ye
     ey_pointv = detadya(j)
     dxy4 = dx1a*dy4a*ey_pointv
     dzy4 = dz1a*dy4a*ey_pointv
     dyy2 = dy2a*ey_pointv
     dyya =  dysa*detady(j+1)*ey_pointv
     dyyam = dysa*detady(j)*ey_pointv
  DO k = zs, ze
  DO i = xs, xe
              
     dqdt(i,k,j,eqn) = dfdt(i,k,j,eqn) + rate(i,k,j,eqn)
    
     dtcpp=dtx(i,k,j)/cp_prndtl
     fxbar = dphidx(i,k)
     fzbar = dphidz(i,k)
     term = one+fxbar**2+fzbar**2
     convx = (uvel(i,k,j))*dx2a+(lambdag(i,k,j+1)-lambdag(i,k,j))*fxbar*dxy4
     convz = (wvel(i,k,j))*dz2a+(lambdag(i,k,j+1)-lambdag(i,k,j))*fzbar*dzy4
     convy = (vvel(i,k,j))*dyy2  &
          + (lbdavg(i,k,j,3)*dphidza(i,k)-lbdavg(i,k-1,j,3)*dphidza(i,k-1))*dzy4  &
          + (lbdavg(i,k,j,2)*dphidxa(i,k)-lbdavg(i-1,k,j,2)*dphidxa(i-1,k))*dxy4

     coeff(2,i,k,j) =  convx + lbdavg(i-1,k,j,2)*dxsa
     coeff(3,i,k,j) = -convx + lbdavg(i,k,j,2) * dxsa

     coeff(4,i,k,j) =  convy + lambdag(i,k,j)*term*dyyam
     coeff(5,i,k,j) = -convy + lambdag(i,k,j+1)*term*dyya

     coeff(6,i,k,j) =  convz + lbdavg(i,k-1,j,3)*dzsa
     coeff(7,i,k,j) = -convz + lbdavg(i,k,j,3)*dzsa

     if(f1st_ord_conv) then   !drop to first order the convective terms
       coeff(2,i,k,j) =  coeff(2,i,k,j) + abs(uvel(i,k,j))*dx2a
       coeff(3,i,k,j) =  coeff(3,i,k,j) + abs(uvel(i,k,j))*dx2a
       coeff(4,i,k,j) =  coeff(4,i,k,j) + abs(vvel(i,k,j))*dyy2
       coeff(5,i,k,j) =  coeff(5,i,k,j) + abs(vvel(i,k,j))*dyy2
       coeff(6,i,k,j) =  coeff(6,i,k,j) + abs(wvel(i,k,j))*dz2a
       coeff(7,i,k,j) =  coeff(7,i,k,j) + abs(wvel(i,k,j))*dz2a
     endif

     if (j == ny-1)  then
        coeff(4,i,k,j)  = lambdag(i,k,j)*term*dyyam + lambdag(i,k,j+1)*term*dyya
        coeff(5,i,k,j)  = zero
     endif

     coeff(2:7,i,k,j) = coeff(2:7,i,k,j)*dtcpp

!..................................................................................     
! TO save memory use Bm to store the orhter terms (corners) in the jacobian Matrix
!..................................................................................     
     coe = dxy4*dtcpp
     Bm(1,1,i,k,j) = -(lambdag(i,k,j)*fxbar + lbdavg(i-1,k,j,2)*dphidxa(i-1,k))*coe
     Bm(1,2,i,k,j) =  (lambdag(i,k,j)*fxbar + lbdavg(i,k,j,2)*dphidxa(i,k))*coe
     Bm(1,3,i,k,j) = -(lambdag(i,k,j+1)*fxbar + lbdavg(i,k,j,2)*dphidxa(i,k))*coe
     Bm(1,4,i,k,j) = - sum(Bm(1,1:3,i,k,j))

     coe = dzy4*dtcpp
     Bm(2,1,i,k,j) = -(lambdag(i,k,j)*fzbar + lbdavg(i,k-1,j,3)*dphidza(i,k-1))*coe
     Bm(2,2,i,k,j) =  (lambdag(i,k,j)*fzbar + lbdavg(i,k,j,3)*dphidza(i,k))*coe
     Bm(2,3,i,k,j) = -(lambdag(i,k,j+1)*fzbar + lbdavg(i,k,j,3)*dphidza(i,k))*coe
     Bm(2,4,i,k,j) = - sum(Bm(2,1:3,i,k,j))

  ENDDO
  ENDDO
  ENDDO   

  endif  !if (eqn.eq.1)

  if (is2D) coeff(6:7,:,:,:) = zero

  j=eqn  !gross
  CALL DIRICHLET_BC(j, finest_mesh)   !gross
  eqn=j  !gross 
  
  DO j = ys, ye
  DO k = zs, ze
  DO i = xs, xe
     coeff(1,i,k,j) = one + sum(coeff(2:7,i,k,j))
     if(j == 1 .and. eqn == 1) then
       coeff(4,i,k,j) = zero
!
       c_elim = Bm(1,1,i,k,j)
       Bm(1,4,i,k,j) = Bm(1,4,i,k,j) + cj2wall*c_elim
       coeff(2,i,k,j) = coeff(2,i,k,j) - cj1wall*c_elim
       dqdt(i,k,j,eqn) = dqdt(i,k,j,eqn) + cjswall*c_elim*uvel(i-1,k,0)
       Bm(1,1,i,k,j) = zero

       c_elim = Bm(1,2,i,k,j)
       Bm(1,3,i,k,j) = Bm(1,3,i,k,j) + cj2wall*c_elim
       coeff(3,i,k,j) = coeff(3,i,k,j) - cj1wall*c_elim
       dqdt(i,k,j,eqn) = dqdt(i,k,j,eqn) + cjswall*c_elim*uvel(i+1,k,0)
       Bm(1,2,i,k,j) = zero

       c_elim = Bm(2,1,i,k,j)
       Bm(2,4,i,k,j) = Bm(2,4,i,k,j) + cj2wall*c_elim
       coeff(6,i,k,j) = coeff(6,i,k,j) - cj1wall*c_elim
       dqdt(i,k,j,eqn) = dqdt(i,k,j,eqn) + cjswall*c_elim*uvel(i,k-1,0)
       Bm(2,1,i,k,j) = zero

       c_elim = Bm(2,2,i,k,j)
       Bm(2,3,i,k,j) = Bm(2,3,i,k,j) + cj2wall*c_elim
       coeff(7,i,k,j) = coeff(7,i,k,j) - cj1wall*c_elim
       dqdt(i,k,j,eqn) = dqdt(i,k,j,eqn) + cjswall*c_elim*uvel(i,k+1,0)
       Bm(2,2,i,k,j) = zero
!       write(*,*)i,k,j, coeff(1,i,k,j), coeff(5,i,k,j)
     endif
  ENDDO
  ENDDO
  ENDDO
 
!--------------------------------------------------------------------
   RETURN
END SUBROUTINE COEFFMATRIXQ
!*********************************************************************


!*********************************************************************
SUBROUTINE RHSQ(eqn)

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

  include "parallel.h"

!---------------------------------------------------------------------
! Dummy variables
  INTEGER, INTENT(IN) :: eqn
! Local variables
! c_i::center point jacobian, b_i and d_i::corners
  INTEGER ::  i, j, k, jp1, js, jmax 
  INTEGER ::  xs,xe,zs,ze,ys,ye,indx,color, col(7)
  REAL*8  ::  coe,c2,c3,c6,c7,b1,b2,b3,b4,d1,d2,d3,d4
  REAL*8  ::  f2,f3,f6,f7,s1,s2,s3,s4,t1,t2,t3,t4
  REAL*8  ::  oldq,rhsg,resid,omg_uvw
!---------------------------------------------------------------------
! Initialize variables
  ys = drange(5)+1
  ye = nyv(eqn)

!  col = (/ 0, 0, eqn, eqn, drange(5), drange(6), 0 /)
  col(3:4)=eqn
  col(5:6) = drange(5:6)
  col(7) = 1
  omg_uvw = 1.00d0

  resid= 0.d0

  jmax = ye+1
  js = ys-1

  do color = 1,2
    col(1:2) = color
    do indx = 1 , eqn_nclr(color,eqn)

       i =  eqn_iclr(indx,color,eqn)
       k =  eqn_kclr(indx,color,eqn)

!--- setup thomas algorithm
       fy(js) = q(i,k,0,eqn)
       ay(js) = 0.0d0
       cy(js) = 0.0d0
       by(js) = 1.0d0
       do j = ys, ye
          jp1 = merge(j+1,j-1, j<ye .OR. eqn<3)
           
          ay(j) = - coeff(4,i,k,j)
          cy(j) = - coeff(5,i,k,j)
          by(j) = coeff(1,i,k,j)

          c2=coeff(2,i,k,j)
          c3=coeff(3,i,k,j)

          b1=Bm(1,1,i,k,j)
          b2=Bm(1,2,i,k,j)
          b3=Bm(1,3,i,k,j)
          b4=Bm(1,4,i,k,j)

          d1=Bm(2,1,i,k,j)
          d2=Bm(2,2,i,k,j)
          d3=Bm(2,3,i,k,j)
          d4=Bm(2,4,i,k,j)
 

          fy(j) =  q(i-1,k,j,eqn)*c2 + q(i+1,k,j,eqn)*c3   &
               + q(i-1,k,j-1,eqn)*b1+q(i+1,k,j-1,eqn)*b2   &
               + q(i+1,k,jp1,eqn)*b3+q(i-1,k,jp1,eqn)*b4   &
               + q(i,k-1,j-1,eqn)*d1+q(i,k+1,j-1,eqn)*d2   &
               + q(i,k+1,jp1,eqn)*d3+q(i,k-1,jp1,eqn)*d4   &
               + dqdt(i,k,j,eqn)

          oldq = fy(j) - ay(j)*q(i,k,j-1,eqn) - cy(j)*q(i,k,j+1,eqn) - by(j)*q(i,k,j,eqn)
          resid = max(resid,abs(oldq))

       end do

       if(eqn < 3) then
          fy(jmax) = q(i,k,nyv(eqn)+1,eqn)
          ay(jmax) = 0.0d0
          cy(jmax) = 0.0d0
          by(jmax) = 1.0d0
       else
          fy(jmax) = 0.0d0
          ay(jmax) = -1.0d0
          cy(jmax) = 0.0d0
          by(jmax) = 1.0d0
       endif

! --- THOMAS algorithm
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
       enddo

    end do  !indx

    CALL parallel_swap4(col,finest_mesh)

  enddo     !color

  if (eqn == 3) then
     do i= drange(1)-1, drange(2)+1
     do k= drange(3)-1, drange(4)+1
        q(i,k,ny,3) = q(i,k,ny-2,3)
     enddo
     enddo
  endif

  call MPI_ALLREDUCE(resid,res_gas,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)

!---------------------------------------------------------------------
  RETURN
END SUBROUTINE RHSQ
!*********************************************************************

