!     **********************************************************************
!
!IMPLICIT BOUNDARY CONDITION SOLVER
!
!     **********************************************************************

 SUBROUTINE BC_IMP(color)
   
   USE GLOBAL_DATA
   USE BC_DATA
   USE MYMPI
   
   IMPLICIT NONE
   
!---------------------------------------------------------------------------
!  Dummy variables
   INTEGER, INTENT(IN) :: color
!  Local variables
   INTEGER :: i,k,iter,itmax,xs,ys,zs,xe,ye,ze,indx,eqn
   REAL*8  :: dlg, const1, const2
   REAL*8  :: rmax,allrmax,prec_gss,ffold
   REAL*8  :: fx, fz,hx,hz,coe_lam1,coe_lam2 
   REAL*8  :: biga,bigb,bigc,bigr,term
!---------------------------------------------------------------------------
   
   xs = drange(1)
   xe = drange(2)
   ys = drange(5)
   ye = drange(6)
   zs = drange(3)
   ze = drange(4)

!...........................................................................
!  ZERO GRADIENT AT NORTH
!...........................................................................
   IF (color ==2 ) THEN
     do k = zs-1, ze+1
     do i = xs-1, xe+1
     do eqn = 1,neqmax-1
        f(i,k,ny,eqn) = f(i,k,ny-1,eqn)
     end do
     end do
     end do
   ENDIF

!
!  TEMPERATURE BOUNDARY CONDITIONS
!
    coe_lam1 = coe_lamb1
    coe_lam2 = coe_lamb2
    eqn =1
    DO indx = 1 , finest_mesh%nclr(color)
       i =  finest_mesh%iclr(indx,color)
       k =  finest_mesh%kclr(indx,color)
       term = capr(1,i,k) - capd(i,k)*f(i,k,1,1) - capdd(i,k)*f(i,k,2,1)  &
            -  cape(i,k)*f(i,k,1,neqmax) - capee(i,k)*f(i,k,2,neqmax)
       fx = f(i+1,k,0,eqn)-f(i-1,k,0,eqn)
       fz = f(i,k+1,0,eqn)-f(i,k-1,0,eqn)
       f(i,k,0,eqn) = (capa(1,i,k)*fx + capb(1,i,k)*fz + term)/capc(1,i,k)
       f(i,k,0,neqmax) = f(i,k,0,eqn)
       lambdag(i,k,0) = coe_lam1*f(i,k,0,1)+coe_lam2
    ENDDO


!
!  SPECIES BOUNDARY CONDITIONS
! 
    hx = 1.0d0/(2.0d0*dx)
    hz = 1.0d0/(2.0d0*dz)

    DO EQN=2,neqmax-1

      DO indx = 1 , finest_mesh%nclr(color)
         i =  finest_mesh%iclr(indx,color)
         k =  finest_mesh%kclr(indx,color)
         term = capr(eqn,i,k) - capd(i,k)*f(i,k,1,eqn) - capdd(i,k)*f(i,k,2,eqn)
         fx=f(i+1,k,0,eqn)-f(i-1,k,0,eqn)
         fz=f(i,k+1,0,eqn)-f(i,k-1,0,eqn)
         f(i,k,0,eqn) = (capa(2,i,k)*fx + capb(2,i,k)*fz + term)/capc(2,i,k)
      ENDDO
    ENDDO
   
!----------------------------------------------------------------------
   RETURN
 END SUBROUTINE BC_IMP
!**********************************************************************
      

!****************************************************************************
!
!THIS SUBROUTINE EVALUATES THE JACOBIANS FOR THE IMPLICIT BOUNDARY CONDITIONS
!
!****************************************************************************
 
 SUBROUTINE BC_JAC

   USE GLOBAL_DATA
   USE BC_DATA
   USE MYMPI

   IMPLICIT NONE

!---------------------------------------------------------------------------
!  Local Variables
   INTEGER ::  i,j,k,eqn, xs,xe,zs,ze
   INTEGER ::  inwt,igss,ngss,nnwt
   INTEGER ::  ys,ye
   REAL*8  :: const4,term
   REAL*8  :: d_dlg,d_rb0,d_vvel0,d_capr,d_capa,d_capb,d_capc
   REAL*8  :: dls,dlg,d_termf
   REAL*8  :: rmax,allrmax
   REAL*8  :: rb0,vvel0
   REAL*8  :: fx, fz
   REAL*8  :: hx,hz
   REAL*8  :: const_der_gas, const_der_sld
!---------------------------------------------------------------------------

!.....................................................................
!  UPDATE THE SURFACE MASS FLUX (PYROLISYS)
!  ONLY IF VELOCITY IS UPDATED. THIS IS NO PROBLEM IN THE OSEEN
!  BUT IN THE NS IT IS DIFFICULT TO UPDATE THE VELOCITY
!.....................................................................
!   CALL PYROLYSIS
!
   hx = 1.0d0/(2.0d0*dx)
   hz = 1.0d0/(2.0d0*dz)
   ys = drange(5)
   ye = drange(6)
   xs = drange(1)
   xe = drange(2)
   zs = drange(3)
   ze = drange(4)
   
      do k=zs,ze
         do i=xs,xe
            term=1.0d0+dphidx(i,k)**2+dphidz(i,k)**2
            dls = lambdas(i,k,0)
            dlg = coe_lamb1*f(i,k,0,1)+coe_lamb2
            d_dlg = coe_lamb1
            rb0 = da_s(i,k,0)*exp(-theta_s(i,k,0)/f(i,k,0,1))
            d_rb0 = rb0*theta_s(i,k,0)/f(i,k,0,1)**2
            vvel0 = rhos(i,k,0)*rb0*sqrt(term)
            vvel(i,k,0) = vvel0
            d_vvel0 = vvel0/rb0*d_rb0
            const4 = term*detady(0)/(2.0d0*dy)

            const_der_gas = const4*dlg
            const_der_sld = const4*dls
 
            capa(1,i,k) = dphidx(i,k)*(dls-dlg)*hx
            capa(2,i,k) = - dphidx(i,k) * dlg * hx
            d_capa = - dphidx(i,k)*d_dlg*hx
            capb(1,i,k) = dphidz(i,k)*(dls-dlg)*hz
            capb(2,i,k) = - dphidz(i,k) * dlg * hz
            d_capb = - dphidz(i,k)*d_dlg*hz

            capc(1,i,k) =  3.0d0*(const_der_gas + const_der_sld)
            d_capc = const4*d_dlg*3.0d0  
            capd(i,k) = -4.0d0*const_der_gas
            capdd(i,k) = +1.0d0*const_der_gas
            cape(i,k) = -4.0d0*const_der_sld
            capee(i,k) = +1.0d0*const_der_sld
            
            capr(1,i,k) = qheats(i,k,0)*vvel0
            d_capr = qheats(i,k,0)*d_vvel0&
                 + const4*d_dlg*(4.0d0*f(i,k,1,1)-f(i,k,2,1))
            fx= (f(i+1,k,0,1)-f(i-1,k,0,1))
            fz= (f(i,k+1,0,1)-f(i,k-1,0,1))
  
            d_termf = -d_capc*f(i,k,0,1)+d_capr+d_capa*fx+d_capb*fz
            
            capr(1,i,k) = capr(1,i,k) - d_termf*f(i,k,0,1)   !SOLVE C*f = A*Fx+B*Fz+R
            do eqn = 2,neqgas
               capr(eqn,i,k) = cp*xsurf(i,k,eqn)*vvel(i,k,0)
            enddo

            capc(1,i,k) = capc(1,i,k) - d_termf              ! C
            capc(2,i,k) = cp*vvel(i,k,0) + 3.0d0*const_der_gas
         enddo
      enddo    
!-------------------------------------------------------------------------------
         RETURN
END SUBROUTINE BC_JAC
!*******************************************************************************
