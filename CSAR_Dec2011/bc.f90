! Heterogeneous Propellant 3D Combustion Code
! Constant Mass Flux
! DOE-ASCI Center for Simulation of Advanced Rockets
! University of Illinois at Urbana-Champaign
! Author: Thomas L. Jackson
! Date last Modified: 16 January 2001
! Parallelized by M. Campbell on May 9, 2001
! //////////////////////////////////////////////////
!
! ********************************************************************
! PROPRIETARY CODE - DO NOT USE WITHOUT PERMISSION OF CSAR!
! ********************************************************************
!
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
!
!     Boundary Conditions
!
!     This function sets the boundary conditions on
!       all four sides of the rectangle.
!
! ********************************************************************
SUBROUTINE BC

USE GLOBAL_DATA
USE BC_DATA
   
IMPLICIT NONE
   
! Local variables
INTEGER :: i, k, xs,zs,xe,ze,eqn
!-------------------------------------------

finest_mesh%w => f
var_updated = 0
   
xs = drange(1)
xe = drange(2)
zs = drange(3)
ze = drange(4)
   
! Grid Boundary Conditions
! North
do eqn = 1,neqmax-1
    do k = zs-1, ze+1
        do i = xs-1, xe+1
            f(i,k,ny,eqn) = f(i,k,ny-1,eqn)
        end do
    end do
end do
! South
do k = zs-1, ze+1
    do i = xs-1, xe+1
        f(i,k,ny,neqmax) = tcold
    end do
end do

CALL BC_TEMP
CALL PYROLYSIS 
CALL BC_SPEC

!----------------------------------------------------------
RETURN
END SUBROUTINE BC
!**************************************************************
      
      
!**************************************************************************
!
!     Solver for species bc  "NEW VERSION"
!
!*************************************************************************

SUBROUTINE BC_SPEC
   
USE GLOBAL_DATA
USE BC_DATA
USE MYMPI
   
IMPLICIT NONE
include "parallel.h"
   
!-----------------------------------------------------------------------------
!  Local variables
INTEGER :: i,k,iter,itmax,xs,zs,xe,ze,eqn
INTEGER :: indx,color,col(7)
REAL*8  :: dlg, const1, const2
REAL*8  :: rmax,allrmax,prec_gss,ffold, addeqn
REAL*8  :: fx, fz,hx,hz,  coe_lam1,coe_lam2
!-----------------------------------------------------------------------------
   
xs = drange(1)
xe = drange(2)
zs = drange(3)
ze = drange(4)
col(3) = 2; col(4) = neqmax-1
col(5) = 0; col(6) = 0
col(7) = 0

hx = 1.0d0/(2.0d0*dx)
hz = 1.0d0/(2.0d0*dz)
coe_lam1 = coe_lamb1
coe_lam2 = coe_lamb2
prec_gss=1.0d-08

do k=zs,ze
    do i=xs,xe
        dlg = coe_lam1*f(i,k,0,1)+coe_lam2
!redefine vvel and lambdag
        vvel(i,k,0) = rhos(i,k,0)*rb(i,k)&
            *sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)
        lambdag(i,k,0) = dlg  
             
        const1 = dlg*(one+dphidx(i,k)**2+dphidz(i,k)**2)&
            *detady(0)/(2.0d0*dy*cp)
        const2 = vvel(i,k,0) + 3.0d0*const1
        capa(2,i,k) = -dlg*dphidx(i,k)/cp*hx
        capc(2,i,k) = one/const2
!         capr(i,k) = const1
        do eqn = 2,neqgas
            capr(eqn,i,k) = xsurf(i,k,eqn)*vvel(i,k,0)+ const1*(4.0d0*f(i,k,1,eqn)-f(i,k,2,eqn))
        enddo
    enddo
enddo
        
itmax = 15
do iter=1,itmax      
    rmax=0.0d0
    DO color = 1,2
        col(1:2) = color
        do indx = 1 , eqn_nclr(color,0)
            i = eqn_iclr(indx,color,0)
            k = eqn_kclr(indx,color,0)
            do eqn=2,neqgas
                fx=f(i+1,k,0,eqn)-f(i-1,k,0,eqn)
                ffold = f(i,k,0,eqn)
                f(i,k,0,eqn) = capc(2,i,k)*(capr(eqn,i,k)+ capa(2,i,k)*fx)
                rmax = max(rmax,abs(f(i,k,0,eqn)-ffold))
            enddo                      !end equations
        enddo                        !end points
        call parallel_swap4(col,finest_mesh)
!       call DIRICHLET_BC(0,finest_mesh)
    ENDDO                          !end colors

    call MPI_ALLREDUCE(rmax,allrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
    if (allrmax.lt.prec_gss) exit
enddo                            !end iterations


!--------------------------------------------------------------W
!        write(*,*)'SPECIES ITER',iter,allrmax,eqn
!--------------------------------------------------------------W
   
if (iter.ge.itmax) then
    write(*,*) 'ERROR @ BC >SPECIES< ITERATION COUNT REACHED'
    write(6,*) iter,rmax
    STOP
endif
   
if(isparallel(3)) then
    if(mywid.eq.0) then
        write(*,*)'Rocfire3D :: Fatal Error!  Please check your yrange!'
    endif
    STOP
endif
   
RETURN
END SUBROUTINE BC_SPEC
!***************************************************************************
      

!**********************************************************************
!
!  Solver for temperature bc  'NEW VERSION'
!
!**********************************************************************
 
SUBROUTINE BC_TEMP

USE GLOBAL_DATA
USE BC_DATA
USE MYMPI

IMPLICIT NONE
include "parallel.h"

!----------------------------------------------------------------------------
!  Local Variables
INTEGER ::  i,k,iter,itmax,xs,xe,zs,ze
INTEGER ::  inwt,igss,ngss,nnwt
INTEGER ::  col(7),indx,color
REAL*8  :: const3,const4,term,termf
REAL*8  :: d_dlg,d_rb0,d_vvel0,d_capr,d_capa,d_capb,d_capc
REAL*8  :: dls,dlg
REAL*8  :: rmax,allrmax
REAL*8  :: rb0,vvel0
REAL*8  :: fx, fz
REAL*8  :: hx,hz
REAL*8  :: coe_lam1,coe_lam2
REAL*8  :: prec_nwt,prec_gss,res_nwt,res_gss,ffold,resid
!----------------------------------------------------------------------------
!
!
hx = 1.0d0/(2.0d0*dx)
hz = 1.0d0/(2.0d0*dz)
xs = drange(1)
xe = drange(2)
zs = drange(3)
ze = drange(4)
col(3) = 1; col(4) = 1
col(5) = 0; col(6) = 0
col(7) = 0
nnwt = 12
ngss = 10
prec_nwt=1.d-08
prec_gss=1.d-08
inwt = 0
res_nwt = 10.0d0
coe_lam1 = coe_lamb1
coe_lam2 = coe_lamb2

!------------------START THE SOLUTION ALGORITHM---------------
NEWTON: DO
    inwt=inwt+1
    resid = zero

    do k=zs,ze
        do i=xs,xe
            term=1.0d0+dphidx(i,k)**2+dphidz(i,k)**2
            dls = lambdas(i,k,0)
            dlg = coe_lam1*f(i,k,0,1)+coe_lam2
            d_dlg = coe_lam1
            rb0 = da_s(i,k,0)*exp(-theta_s(i,k,0)/f(i,k,0,1))
            d_rb0 = rb0*theta_s(i,k,0)/f(i,k,0,1)**2
            vvel0 = rhos(i,k,0)*rb0*sqrt(term)
            d_vvel0 = vvel0/rb0*d_rb0
            const4 = term*detady(0)/(2.0d0*dy)
            const3 = const4*(dlg+dls)*3.0d0
            d_capc = const4*d_dlg*3.0d0
            capa(1,i,k) = dphidx(i,k)*(dls-dlg)*hx
            d_capa = - dphidx(i,k)*d_dlg*hx
            capb(1,i,k) = dphidz(i,k)*(dls-dlg)*hz
            d_capb = - dphidz(i,k)*d_dlg*hz
            capc(1,i,k) = const3
            capr(1,i,k) = qheats(i,k,0)*vvel0&
                 + const4*dlg*(4.0d0*f(i,k,1,1)-f(i,k,2,1))&
                 + const4*dls*(4.0d0*f(i,k,1,neqmax)-f(i,k,2,neqmax))
            d_capr = qheats(i,k,0)*d_vvel0&
                 + const4*d_dlg*(4.0d0*f(i,k,1,1)-f(i,k,2,1))
            fx=(f(i+1,k,0,1)-f(i-1,k,0,1))
            fz=(f(i,k+1,0,1)-f(i,k-1,0,1))
            resid = max(resid,abs(capa(1,i,k)*fx+capb(1,i,k)*fz&
                    -capc(1,i,k)*f(i,k,0,1)+capr(1,i,k)))
            termf = -d_capc*f(i,k,0,1)+d_capr+d_capa*fx+d_capb*fz
            capr(1,i,k) = capr(1,i,k) - termf*f(i,k,0,1)   !SOLVE C*f = A*Fx+B*Fz+R
            capc(1,i,k) = 1.0d0/(capc(1,i,k) - termf )  ! 1/C
        enddo
    enddo

    call MPI_ALLREDUCE(resid,res_nwt,1,MPI_DOUBLE_PRECISION,&
        MPI_MAX,comm3d,ierr)
    IF(res_nwt <= prec_nwt .OR. inwt >= nnwt) EXIT NEWTON

    igss=0
    GAUSS: DO
        igss = igss+1
        rmax = 0
         
        DO color = 1,2
            col(1:2) = color
            DO indx = 1 , eqn_nclr(color,0)
                i = eqn_iclr(indx,color,0)
                k = eqn_kclr(indx,color,0)
                fx=f(i+1,k,0,1)-f(i-1,k,0,1)
                fz=f(i,k+1,0,1)-f(i,k-1,0,1)
                ffold = f(i,k,0,1)
                f(i,k,0,1) = capc(1,i,k)*(capr(1,i,k)+capa(1,i,k)*fx+capb(1,i,k)*fz)
                rmax=rmax+abs(f(i,k,0,1)-ffold)
            ENDDO
            call parallel_swap4(col,finest_mesh)
!           call DIRICHLET_BC(0,finest_mesh)
        ENDDO

        call MPI_ALLREDUCE(rmax,res_gss,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
        if(res_gss <= prec_gss .OR. igss >= ngss) EXIT GAUSS

    END DO GAUSS
END DO NEWTON
   
!
!  IMPOSE the continuity of the temperature across the interface
!
f(:,:,0,neqmax) = f(:,:,0,1)

if(inwt.gt.nnwt) then
    write(6,*) 'Rocfire3D :: ERROR BC_TEMP @ NEWTON',inwt,res_nwt,igss
endif

if(isparallel(3)) then
    write(*,*)'Rocfire3D :: You need to Broadcast f(i,k,0,1) to colcom'
    write(*,*)'Rocfire3D :: FATAL ERROR! Goodbye.'
    STOP       
endif
         
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE BC_TEMP
!*******************************************************************************
