
!****************************************************************************
!
!THIS SUBROUTINE EVALUATES THE JACOBIANS FOR THE IMPLICIT BOUNDARY CONDITIONS
!
!****************************************************************************
 
 SUBROUTINE BC_JAC

   USE GLOBAL_DATA
   USE implicit
   USE BC_DATA
   USE MYMPI

   IMPLICIT NONE

!---------------------------------------------------------------------------
!  Local Variables
   INTEGER ::  i,j,k,eqn, xs,xe,zs,ze
   INTEGER ::  inwt,igss,ngss,nnwt
   INTEGER ::  ys,ye
   REAL*8  :: const4,term
   REAL*8  :: d_dlg,d_rb0,d_vvel0,d_capr,d_capa,d_capb,d_capc,d_capd,d_capdd
   REAL*8  :: dls,dlg,d_termf
   REAL*8  :: rmax,allrmax
   REAL*8  :: rb0,vvel0
   REAL*8  :: fx, fz
   REAL*8  :: hx,hz
   REAL*8  :: const_der_gas, const_der_sld,d_const_der_gas
   REAL*8  :: xsnew(4),dxsnew(4)
   REAL*8  :: tmp_f,d_tmp_f,qheat_p1,d_qheat_p1,qheat0,d_qheat0,tsinv2
!---------------------------------------------------------------------------

!.....................................................................
!  UPDATE THE SURFACE MASS FLUX (PYROLISYS)
!  ONLY IF VELOCITY IS UPDATED. THIS IS NO PROBLEM IN THE OSEEN
!  BUT IN THE NS IT IS DIFFICULT TO UPDATE THE VELOCITY
!.....................................................................
!   CALL PYROLYSIS
!
!
!  TLJ
   if(.not. implicitBC)then  !cancel <MOVEUP>
      dfdt(:,:,0,1:maxsize) = 0d0
      return
   endif   !cancel <MOVEUP>
!
   ys = drange(5)
   ye = drange(6)
   xs = drange(1)
   xe = drange(2)
   zs = drange(3)
   ze = drange(4)

   hx = 1.0d0/(2.0d0*dx)
   hz = 1.0d0/(2.0d0*dz)
   kloop: do k=zs,ze
      iloop: do i=xs,xe

         term=1.0d0+dphidx(i,k)**2+dphidz(i,k)**2
         dls = lambdas(i,k,0)

         dlg = a_lamg*f(i,k,0,1)**e_lamg + b_lamg
         d_dlg = e_lamg*a_lamg*f(i,k,0,1)**(e_lamg-one)

         rb0 = da_s(i,k,0)*exp(-theta_s(i,k,0)/f(i,k,0,1))
         d_rb0 = rb0*theta_s(i,k,0)/f(i,k,0,1)**2

         vvel0 = rhos(i,k,0)*rb0*sqrt(term)
         vvel(i,k,0) = vvel0
         d_vvel0 = vvel0/rb0*d_rb0

         const4 = term*detady(0)/(2.0d0*dy)

         tsinv2 = one/(f(i,k,0,1)*f(i,k,0,1))
         tmp_f = f_pyr(1)*exp(-f_pyr(2)/f(i,k,0,1))*press**f_pyr(3)*f(i,k,0,1)**f_pyr(4)
         tmp_f = tmp_f*psurf(i,k,2)
         d_tmp_f = f_pyr(2)*tmp_f*tsinv2 + f_pyr(4)*tmp_f/f(i,k,0,1)
         qheat_p1 = qheat_OX + tmp_f *Qg(1)
         d_qheat_p1 = d_tmp_f *Qg(1)
         qheat0 = psurf(i,k,3)*qheat_p1 + (1.0d0-psurf(i,k,3))*qheat_binder
         d_qheat0 = psurf(i,k,3)*d_qheat_p1
!
         xsnew(2) = psurf(i,k,3)*(1d0 - tmp_f)
         xsnew(3) = psurf(i,k,3)*tmp_f
         xsnew(4) = 1.0d0-psurf(i,k,3)
         dxsnew(2) = psurf(i,k,3)*(- d_tmp_f)
         dxsnew(3) = psurf(i,k,3)*d_tmp_f
         dxsnew(4) = zero     

!        TLJ
         qheat0 = qheats(i,k,0)
         d_qheat0 = zero
         xsnew(2) = xsurf(i,k,2)
         xsnew(3) = xsurf(i,k,3)
         xsnew(4) = xsurf(i,k,4)
         dxsnew = zero
!
         const_der_gas = const4*dlg
         const_der_sld = const4*dls
         d_const_der_gas = const4*d_dlg

         capa(1,i,k) = dphidx(i,k)*(dls-dlg)*hx
         capa(2,i,k) = - dphidx(i,k) * dlg * hx
         d_capa = - dphidx(i,k)*d_dlg*hx

         capb(1,i,k) = dphidz(i,k)*(dls-dlg)*hz
         capb(2,i,k) = - dphidz(i,k) * dlg * hz
         d_capb = - dphidz(i,k)*d_dlg*hz

         capc(1,i,k) =  3.0d0*(const_der_gas + const_der_sld)
         capc(2,i,k) =  3.0d0*const_der_gas     !discr  A*Fx+B*Fz+-C*F+R
         d_capc = const4*d_dlg*3.0d0  

         capd(i,k) = -4.0d0*const_der_gas
         capdd(i,k) = +1.0d0*const_der_gas
         d_capd = -4.0d0*d_const_der_gas
         d_capdd = +1.0d0*d_const_der_gas

         cape(i,k) = -4.0d0*const_der_sld
         capee(i,k) = +1.0d0*const_der_sld

         fx=(f(i+1,k,0,1)-f(i-1,k,0,1))
         fz=(f(i,k+1,0,1)-f(i,k-1,0,1))
         dfdt(i,k,0,1) = qheat0*vvel0 &   !SOLVE C*f = A*Fx+B*Fz+R
              &- capd(i,k)*f(i,k,1,1) - capdd(i,k)*f(i,k,2,1) &
              &- cape(i,k)*f(i,k,1,neqmax) - capee(i,k)*f(i,k,2,neqmax) &
              &+ capa(1,i,k)*fx + capb(1,i,k)*fz - capc(1,i,k) * f(i,k,0,1)

         d_capr = qheat0*d_vvel0 + vvel0*d_qheat0 &
              & - d_capd*f(i,k,1,1) - d_capdd*f(i,k,2,1) &
              & - d_capc*f(i,k,0,1) + d_capa*fx + d_capb*fz

         capc(1,i,k) = capc(1,i,k) - d_capr        ! C
         capc(2,i,k) = capc(2,i,k) + cp*vvel0
         d_capc = d_capc + cp*d_vvel0

!species temperature correlation
         do eqn = 2,4
            fx= f(i+1,k,0,eqn)-f(i-1,k,0,eqn)
            fz= f(i,k+1,0,eqn)-f(i,k-1,0,eqn)
            capT(eqn,i,k) = cp*vvel0*dxsnew(eqn) +  cp*xsnew(eqn)*d_vvel0 &
                 &- d_capd*f(i,k,1,eqn) - d_capdd*f(i,k,2,eqn) & 
                 &+ d_capa*fx + d_capb*fz - d_capc* f(i,k,0,eqn) 
            dfdt(i,k,0,eqn)  = cp*xsnew(eqn)*vvel(i,k,0) &   !SOLVE C*f = A*Fx+B*Fz+R
                 &- capd(i,k)*f(i,k,1,eqn) - capdd(i,k)*f(i,k,2,eqn) & 
                 &+ capa(2,i,k)*fx + capb(2,i,k)*fz - capc(2,i,k) * f(i,k,0,eqn)
         enddo
      enddo iloop
   enddo kloop
   dfdt(:,:,0,maxsize) = 0d0
!--------------------- ---------------------------------------------------------
   RETURN
 END SUBROUTINE BC_JAC
!*******************************************************************************
