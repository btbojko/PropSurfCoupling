! ********************************************************************
! ********************************************************************
! Updated: TLJ; 12/19/2016
! Filename: bcnew.f90
! ********************************************************************
! ********************************************************************
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
!
!     This function sets the boundary conditions for T_gas,
!     T_solid, and Y_gas,i on the North, South, and Interface.
!     At the interface it uses Red-Black Gauss-Seidel to
!     solve the linear system.
!
!     subroutines: BC
!                  BC_SPEC
!                  BC_TEMP
!                  BC_TCHECK
!
!
! ********************************************************************
SUBROUTINE BC

  USE GLOBAL_DATA
  USE BC_DATA

  IMPLICIT NONE

! Local variables
  INTEGER :: i, k, xs,zs,xe,ze,eqn

  finest_mesh%w => f
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
! Interface
  CALL BC_TEMP
  CALL SOLID_PROP
  call PYROLYSIS
  CALL BC_SPEC
  if (NORADIATION.eqv..FALSE.) call BC_TCHECK

  RETURN
END SUBROUTINE BC


!*********************************************************************
!
!     Solver for species bc  "NEW VERSION"
!
!*********************************************************************

SUBROUTINE BC_SPEC

  USE GLOBAL_DATA
  USE BC_DATA
  USE MYMPI

  IMPLICIT NONE
  include "parallel.h"

!---------------------------------------------------------------------
!  Local variables
  INTEGER :: i,k,iter,itmax,xs,zs,xe,ze,eqn
  INTEGER :: indx,color,col(7)
  REAL*8  :: dlg, const1, const2
  REAL*8  :: rmax,allrmax,prec_gss,ffold, addeqn
  REAL*8  :: fx, fz,hx,hz
!---------------------------------------------------------------------

  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  col(3) = 2
  col(4) = neqmax-1 ! Needed for Red-Black Gauss-Seidel Solver
  col(5) = 0
  col(6) = 0
  col(7) = 1

  hx = 1.0d0/(2.0d0*dx)
  hz = 1.0d0/(2.0d0*dz)
  prec_gss=1.0d-08

  do k=zs,ze
  do i=xs,xe
  do eqn = 2,neqgas
     dlg = lambdag(i,k,0)/lewis(eqn)
     const1 = dlg*(one+dphidx(i,k)**2+dphidz(i,k)**2)&
                   *detady(0)/(2.0d0*dy*cp)
     const2 = vvel(i,k,0) + 3.0d0*const1

     capa(eqn,i,k) = -dlg*dphidx(i,k)/cp*hx
     capb(eqn,i,k) = -dlg*dphidz(i,k)/cp*hz
     capc(eqn,i,k) = one/const2
     capr(eqn,i,k) = xsurf(i,k,eqn)*vvel(i,k,0)&
                + const1*(4.0d0*f(i,k,1,eqn)-f(i,k,2,eqn))
  enddo
  enddo
  enddo


  itmax = 15
  do iter=1,itmax
     rmax=0.0d0
     DO color = 1,2
        col(1:2) = color
        do indx = 1 , finest_mesh%nclr(color)
           i =  finest_mesh%iclr(indx,color)
           k =  finest_mesh%kclr(indx,color)
           do eqn=2,neqgas
              fx = f(i+1,k,0,eqn)-f(i-1,k,0,eqn)
              fz = f(i,k+1,0,eqn)-f(i,k-1,0,eqn)
              ffold = f(i,k,0,eqn)
              f(i,k,0,eqn) = capc(eqn,i,k)*(capr(eqn,i,k) &
                   + capa(eqn,i,k)*fx+capb(eqn,i,k)*fz)
              rmax = max(rmax,abs(f(i,k,0,eqn)-ffold))
           enddo                     !end equations
        enddo                        !end points
        call parallel_swap4(col,finest_mesh)
     ENDDO                           !end colors

     call MPI_ALLREDUCE(rmax,allrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
     if (allrmax.lt.prec_gss) exit

  enddo  !end iterations

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
!*********************************************************************


!     ****************************************************************
!
!     Solver for temperature bc  'NEW VERSION'
!
!     ****************************************************************

SUBROUTINE BC_TEMP

  USE GLOBAL_DATA
  USE BC_DATA
  USE MYMPI
  use implicit, only : implicitbc

  IMPLICIT NONE
  include "parallel.h"

!---------------------------------------------------------------------
!  Local Variables
  INTEGER :: i,k,j,n,iter,itmax,xs,xe,zs,ze
  INTEGER :: inwt,igss,ngss,nnwt,mynclip
  INTEGER :: col(7),indx,color,wunit
  REAL*8  :: const3,const4,term,termf
  REAL*8  :: d_dlg,d_rb0,d_vvel0,d_capr,d_capa,d_capb,d_capc
  REAL*8  :: dls,dlg
  REAL*8  :: rmax,allrmax
  REAL*8  :: rb0,vvel0
  REAL*8  :: fx, fz
  REAL*8  :: hx,hz
  REAL*8  :: prec_nwt,prec_gss,res_nwt,res_gss,ffold,resid
  REAL*8  :: tmp_f,d_tmp_f,qheat_p1,d_qheat_p1,qheat0,d_qheat0,tsinv2
  logical :: outcome
!---------------------------------------------------------------------
!
!
  hx = 1.0d0/(2.0d0*dx)
  hz = 1.0d0/(2.0d0*dz)
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  col(3) = 1
  col(4) = 1   ! Needed for Red-Black Gauss-Seidel Solver
  col(5) = 0
  col(6) = 0
  col(7) = 1
  nnwt = 12
  ngss = 10
  prec_nwt=1.d-08
  prec_gss=1.d-08
  inwt = 0
  res_nwt = 10.0d0


!------------------START THE SOLUTION ALGORITHM---------------
  j = 0
  NEWTON: DO
     inwt=inwt+1
     resid = zero

     if (NORADIATION) then
        do k=zs,ze
        do i=xs,xe

           term=1.0d0+dphidx(i,k)**2+dphidz(i,k)**2

           dls = lambdas(i,k,0)
           dlg = a_lamg*f(i,k,0,1)**e_lamg + b_lamg
           d_dlg = e_lamg*a_lamg*f(i,k,0,1)**(e_lamg-one)

           rb0 = da_s(i,k,0)*exp(-theta_s(i,k,0)/f(i,k,0,1))
           d_rb0 = rb0*theta_s(i,k,0)/f(i,k,0,1)**2

           vvel0 = rhos(i,k,0)*rb0*sqrt(term)
           d_vvel0 = vvel0/rb0*d_rb0

           const4 = term*detady(0)/(2.0d0*dy)
           const3 = const4*(dlg+dls)*3.0d0

           capa(1,i,k) = dphidx(i,k)*(dls-dlg)*hx
           d_capa = - dphidx(i,k)*d_dlg*hx

           capb(1,i,k) = dphidz(i,k)*(dls-dlg)*hz
           d_capb = - dphidz(i,k)*d_dlg*hz

           capc(1,i,k) = const3
           d_capc = const4*d_dlg*3.0d0

           capr(1,i,k) = qheats(i,k,0)*vvel0&
                + const4*dlg*(4.0d0*f(i,k,1,1)-f(i,k,2,1))&
                + const4*dls*(4.0d0*f(i,k,1,neqmax)-f(i,k,2,neqmax))
           capr(1,i,k) = capr(1,i,k) + BCradheat(i,k)  !radiation
           d_capr = qheats(i,k,0)*d_vvel0 &
                + const4*d_dlg*(4.0d0*f(i,k,1,1)-f(i,k,2,1))

           fx = (f(i+1,k,0,1)-f(i-1,k,0,1))
           fz = (f(i,k+1,0,1)-f(i,k-1,0,1))

           resid = max(resid,abs(capa(1,i,k)*fx+capb(1,i,k)*fz&
                -capc(1,i,k)*f(i,k,0,1)+capr(1,i,k)))

           termf = -d_capc*f(i,k,0,1)+d_capr+d_capa*fx+d_capb*fz
           capr(1,i,k) = capr(1,i,k) - termf*f(i,k,0,1)   !SOLVE C*f = A*Fx+B*Fz+R
           capc(1,i,k) = 1.0d0/(capc(1,i,k) - termf )  ! 1/C

        enddo
        enddo

     else

        do k=zs,ze
        do i=xs,xe

           term=1.0d0+dphidx(i,k)**2+dphidz(i,k)**2

           dls = lambdas(i,k,0)
           dlg = a_lamg*f(i,k,0,1)**e_lamg + b_lamg
           d_dlg = e_lamg*a_lamg*f(i,k,0,1)**(e_lamg-one)

           rb0 = da_s(i,k,0)*exp(-theta_s(i,k,0)/f(i,k,0,1))
           d_rb0 = rb0*theta_s(i,k,0)/f(i,k,0,1)**2
           vvel0 = rhos(i,k,0)*rb0*sqrt(term)
           d_vvel0 = vvel0/rb0*d_rb0

           tsinv2 = one/(f(i,k,0,1)*f(i,k,0,1))
           tmp_f = f_pyr(1)*exp(-f_pyr(2)/f(i,k,0,1))*press**f_pyr(3)*f(i,k,0,1)**f_pyr(4)
           tmp_f = tmp_f*psurf(i,k,2)
           d_tmp_f = f_pyr(2)*tmp_f*tsinv2 + f_pyr(4)*tmp_f/f(i,k,0,1)
           qheat_p1 = qheat_OX + tmp_f *Qg(1)
           d_qheat_p1 = d_tmp_f *Qg(1)
           qheat0 = (psurf(i,k,3)*qheat_p1 + (1.0d0-psurf(i,k,3))*qheat_binder)
           d_qheat0 = psurf(i,k,3)*d_qheat_p1

           const4 = term*detady(0)/(2.0d0*dy)
           const3 = const4*(dlg+dls)*3.0d0
           d_capc = const4*d_dlg*3.0d0
           capa(1,i,k) = dphidx(i,k)*(dls-dlg)*hx
           d_capa = - dphidx(i,k)*d_dlg*hx
           capb(1,i,k) = dphidz(i,k)*(dls-dlg)*hz
           d_capb = - dphidz(i,k)*d_dlg*hz
           capc(1,i,k) = const3
           capr(1,i,k) = qheat0*vvel0&
                + const4*dlg*(4.0d0*f(i,k,1,1)-f(i,k,2,1))&
                + const4*dls*(4.0d0*f(i,k,1,neqmax)-f(i,k,2,neqmax))
           capr(1,i,k) = capr(1,i,k) + BCradheat(i,k)  !radiation
           d_capr = qheat0*d_vvel0 + vvel0*d_qheat0&
                + const4*d_dlg*(4.0d0*f(i,k,1,1)-f(i,k,2,1))
           fx=(f(i+1,k,0,1)-f(i-1,k,0,1))
           fz=(f(i,k+1,0,1)-f(i,k-1,0,1))
           resid = max(resid,abs(capa(1,i,k)*fx+capb(1,i,k)*fz&
                -capc(1,i,k)*f(i,k,0,1)+capr(1,i,k)))
           termf = -d_capc*f(i,k,0,1)+d_capr+d_capa*fx+d_capb*fz
           capr(1,i,k) = capr(1,i,k) - termf*f(i,k,0,1)   !SOLVE C*f = A*Fx+B*Fz+R
           capc(1,i,k) = 1.0d0/(capc(1,i,k) - termf )  ! 1/C

           !reassign the new evaluated variables
           lambdag(i,k,j) = dlg
           vvel(i,k,j) = vvel0
           qheats(i,k,j) = qheat0

        enddo
        enddo
     endif

     call MPI_ALLREDUCE(resid,res_nwt,1,MPI_DOUBLE_PRECISION,&
          MPI_MAX,comm3d,ierr)
     IF(res_nwt <= prec_nwt .OR. inwt >= nnwt) EXIT NEWTON

     igss=0
     GAUSS: DO
        igss = igss+1
        rmax = 0

        DO color = 1,2
           col(1:2) = color
           DO indx = 1 , finest_mesh%nclr(color)
              i =  finest_mesh%iclr(indx,color)
              k =  finest_mesh%kclr(indx,color)
              fx=f(i+1,k,0,1)-f(i-1,k,0,1)
              fz=f(i,k+1,0,1)-f(i,k-1,0,1)
              ffold = f(i,k,0,1)
              f(i,k,0,1) = capc(1,i,k)*(capr(1,i,k)+capa(1,i,k)*fx+capb(1,i,k)*fz)
              rmax=rmax+abs(f(i,k,0,1)-ffold)
           ENDDO
           call parallel_swap4(col,finest_mesh)
        ENDDO

        call MPI_ALLREDUCE(rmax,res_gss,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
        if(res_gss <= prec_gss .OR. igss >= ngss) EXIT GAUSS

     END DO GAUSS
     IF(inwt >= nnwt) EXIT NEWTON

  END DO NEWTON


!
!  IMPOSE the continuity of the temperature across the interface
!
  f(:,:,0,neqmax) = f(:,:,0,1)


  if(inwt.gt.nnwt) then
     write(6,*) 'Rocfire3D :: ERROR BC_TEMP @ NEWTON',inwt,res_nwt,igss
  endif

  if(isparallel(3)) then
     write(*,*)'Rocfire3D :: FATAL ERROR! Goodbye.'
     STOP       
  endif

!---------------------------------------------------------------------
  RETURN
END SUBROUTINE BC_TEMP
!*********************************************************************
SUBROUTINE BC_TCHECK

  USE GLOBAL_DATA
  USE BC_DATA
  USE MYMPI
  use implicit, only : implicitbc

  IMPLICIT NONE
  include "parallel.h"

!---------------------------------------------------------------------
! Local Variables
  INTEGER ::  i,k,j,n,iter,itmax,xs,xe,zs,ze,mynclip
  real*8 :: ffold
!---------------------------------------------------------------------

  n = 1;j=0
  do k = lbound(f,2),ubound(f,2)
  do i = lbound(f,1),ubound(f,1)
     ffold = f(i,k,j,n)
     f(i,k,0,1) = max( min(f(i,k,0,1),&
             &merge(Tlimit(1),Tlimit(2),psi(i,k,0)>0d0&
             & .or. alp_V(1) > 0.8d0) ),&
             &merge(3d2,3d2,psi(i,k,0)>0d0.or. alp_V(1) > 0.5d0))
     if(f(i,k,j,n) /=  ffold) mynclip = mynclip + 1
     f(i,k,0,neqmax) = f(i,k,0,1)
  enddo
  enddo
!
!
end SUBROUTINE BC_TCHECK
