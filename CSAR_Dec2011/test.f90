! ********************************************************************
SUBROUTINE RHS_test

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!..............................................................................
! Local variables
  INTEGER ::  i, j, jp, k, eqn, ys,ye,xs,xe,zs,ze
  REAL*8  :: dx1, dx2, dy1, dy2, dz1, dz2, dx1a, dx2a, dy1a, dy2a, dz1a, dz2a,term
  REAL*8  :: dxy1, dzy1,dlx1,dlx2,dlz1,dlz2,dly1,dly2,gx,gy,gxx,gyy,gxy,gz,gzz,gzy
  REAL*8  :: conv,diff,myres(2),tmp_res
!..............................................................................

!
! Initialize variables
!
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

!  dfdt = zero

!!CALL BC

!!  f(:,:,ny,1:neqmax-1) = f(:,:,ny-2,1:neqmax-1)

  do eqn=1,neqmax
!   CALL FILL_F_GHOST(eqn)
  enddo

  ys = drange(5)
  ye = drange(6)
  if(ddrange(5).eq.0) ys = ys + 1
  if(ddrange(6).eq.ny) ye = ye - 1
  zs = drange(3)
  ze = drange(4)
  xs = drange(1)
  xe = drange(2)

!  CALL VELOCITY(0)
!  CALL LAMBDA(0)  
  CALL RATE_EXP

  myres = zero
! Loop over Gas-phase equations
  do eqn=1,neqmax-1

     do j = ys, ye
        do k = zs, ze
           do i = xs, xe

              gx = (f(i+1,k,j,eqn)-f(i-1,k,j,eqn))*dx1a
              gy = (f(i,k,j+1,eqn)-f(i,k,j-1,eqn))*dy1a
              gz = (f(i,k+1,j,eqn)-f(i,k-1,j,eqn))*dz1a

              gxx = (f(i+1,k,j,eqn)-2.0d0*f(i,k,j,eqn)&
                   +f(i-1,k,j,eqn))*dx2a
              gyy = (f(i,k,j+1,eqn)-2.0d0*f(i,k,j,eqn)&
                   +f(i,k,j-1,eqn))*dy2a
              gzz = (f(i,k+1,j,eqn)-2.0d0*f(i,k,j,eqn)&
                   +f(i,k-1,j,eqn))*dz2a

              gxy = (f(i+1,k,j+1,eqn)-f(i+1,k,j-1,eqn)&
                   -f(i-1,k,j+1,eqn)+f(i-1,k,j-1,eqn))*dxy1
              gzy = (f(i,k+1,j+1,eqn)-f(i,k+1,j-1,eqn)&
                   -f(i,k-1,j+1,eqn)+f(i,k-1,j-1,eqn))*dzy1

              term=1.0d0+dphidx(i,k)**2+dphidz(i,k)**2

              conv = - vvel(i,k,j)*detady(j)*gy &
                   - uvel(i,k,j)*gx-wvel(i,k,j)*gz
              diff =  lambdag(i,k,j)*(gxx + gzz +&
                   (term)*(&
                   (detady(j)**2)*gyy+deta2dy2(j)*gy)&
                   -2.0d0*dphidx(i,k)*detady(j)*gxy&
                   -2.0d0*dphidz(i,k)*detady(j)*gzy&
                   -(dphi2dx2(i,k)+dphi2dz2(i,k))*detady(j)*gy)/cp&
                   +((dlgdx(i,k,j)-dphidx(i,k)*detady(j)*dlgdy(i,k,j))*&
                   (gx-dphidx(i,k)*detady(j)*gy)&
                   +(dlgdz(i,k,j)-dphidz(i,k)*detady(j)*dlgdy(i,k,j))*&
                   (gz-dphidz(i,k)*detady(j)*gy)&
                   +(detady(j)**2)*dlgdy(i,k,j)*gy)/cp
              tmp_res = (conv+diff+rate(i,k,j,eqn))
!               dfdt(i,k,j,eqn) = (conv+diff+rate(i,k,j,eqn))
              myres(1) = max(myres(1),abs(tmp_res))
              if(i+k+eqn == 1)write(14,*)j,tmp_res,'GAS',conv,diff,rate(i,k,j,eqn),&
                   lambdag(i,k,j),vvel(i,k,j),detady(j)*gy,dy1a,f(i,k,j+1,eqn),f(i,k,j-1,eqn)

           end do
        end do
     end do
  end do

  eqn = neqmax

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

           term=1.0d0+dphidx(i,k)**2+dphidz(i,k)**2

           conv = + (vvel(i,k,0)*detady(j)*gy) &
                *rhos(i,k,j)/rhos(i,k,0)
           diff = (gxx + gzz +&
                (term)*(&
                (detady(j)**2)*gyy+deta2dy2(j)*lambdas(i,k,j)*gy)&
                + dphidx(i,k)*detady(j)*gxy&
                + dphidz(i,k)*detady(j)*gzy&
                + (dphi2dx2(i,k)+dphi2dz2(i,k))*detady(j)*&
                lambdas(i,k,j)*gy)/cp
           tmp_res = conv+diff
           if(i+k == 0)write(14,*)j,conv,diff,rate(i,k,j,eqn),gy,vvel(i,k,0)
           myres(2) = max(myres(2),abs(tmp_res))
        end do
     end do
  end do

  call MPI_ALLREDUCE(myres,maxres,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)

  if (myid == 0) write(*,*)'TEST RESIDUAL',maxres,'NCYC',ncyc,'DT=',dt
  close(14)
!  print*,myid,'SOLID RES',maxval(abs(dfdt(:,:,:,neqmax)))

!--------------------------------------------------------------------------
  RETURN
END SUBROUTINE RHS_test
!*********************************************************************


!**************************************************************************

SUBROUTINE RHSQ_test(eqn)

  USE GLOBAL_DATA
  IMPLICIT NONE

!---------------------------------------------------------------------
! Local variables
  INTEGER ::  i, j, k, eqn, ys,ye,xs,xe,zs,ze, jp
  REAL*8  :: dx1, dx2, dy1, dy2, dz1, dz2, dx1a, dx2a, dy1a, dy2a, dz1a, dz2a,term
  REAL*8  :: dyy,gx,gy,gz, conv, diff, prndtl_cp, maxreq
  REAL*8  :: ey,firstcorr,corrflux,nu,s1,s2,s3,sig1,sig2,dya,gyy
!---------------------------------------------------------------------
  prndtl_cp = pr/cp

! Initialize variables
  dx1 = dx
  dy1 = dy
  dz1 = dz

  dx1a = 1.0d0/dx1
  dx2a = 1.0d0/(2.0d0*dx1)
  dy1a = 1.0d0/dy1
  dy2a = 1.0d0/(2.0d0*dy1)
  dz1a = 1.0d0/dz1
  dz2a = 1.0d0/(2.0d0*dz1)

  ys = drange(5)
  ye = drange(6)
  if(ddrange(5).eq.0) ys = ys + 1
  if(ddrange(6).eq.ny) ye = ye - 1
  zs = drange(3)
  ze = drange(4)
  xs = drange(1)
  xe = drange(2)

!
!  PUT zero here instead of  uvel, wvel beacuse that contribution has already been
!  added to the residual vector
!
  do i=xs-1,xe+1
  do k=zs-1,ze+1
     q(i,k,0,1)=(zero-6.0*q(i,k,1,1)+q(i,k,2,1))*third
     q(i,k,0,2)=(zero-6.0*q(i,k,1,2)+q(i,k,2,2))*third
  enddo
  enddo

  CALL VISCFLUX(eqn)

    maxreq = zero
    do j = ys, ye
     ey = merge(detady(j),half*(detady(j)+detady(j+1)),eqn <3)
     dyy=dy1a*ey
     do k = zs, ze
        do i = xs, xe

           gyy = (q(i,k,j+1,eqn)-two*q(i,k,j,eqn)+q(i,k,j-1,eqn))

           nu = abs(vvel(i,k,j))*dy2a*ey
           firstcorr = -gyy*nu
           jp = j + merge(0,1,vvel(i,k,j)>0)
           s1 = (q(i,k,jp+1,eqn) - q(i,k,jp,eqn))
           s2 = (q(i,k,jp,eqn) - q(i,k,jp-1,eqn))
           s3 = (q(i,k,jp-1,eqn) - q(i,k,jp-2,eqn))
           if(jp==1) s3=s2;
           if(jp==ny) s1=s2

           sig1 = (sign(half,s1)+sign(half,s2))*min(abs(s1),abs(s2))
           sig2 = (sign(half,s2)+sign(half,s3))*min(abs(s2),abs(s3))

           corrflux = nu * (sig1-sig2)

           gx = (q(i+1,k,j,eqn)-q(i-1,k,j,eqn))*dx2a
           gz = (q(i,k+1,j,eqn)-q(i,k-1,j,eqn))*dz2a
           gy = (q(i,k,j+1,eqn)-q(i,k,j-1,eqn))*dy2a

           conv = - vvel(i,k,j)*gy*ey   &            !- firstcorr - corrflux&
                  - uvel(i,k,j)*gx-wvel(i,k,j)*gz
           diff = (fdiff(i,k,j,1)-fdiff(i-1,k,j,1))*dx1a&
                 +(fdiff(i,k,j,2)-fdiff(i,k-1,j,2))*dz1a&
                 +(fdiff(i,k,j,3)-fdiff(i,k,j-1,3))*dyy
!......................................................................
! TO improve efficiency dqdt was previously premultiplied by dtx
!......................................................................
           dqdt(i,k,j,eqn) =  (dqdt(i,k,j,eqn) - q(i,k,j,eqn))/dtx(i,k,j)+conv+prndtl_cp*diff
!!$           dqdt(i,k,j,eqn) = q(i,k,j,eqn)-dtx(i,k,j)*dqdt(i,k,j,eqn) 
           maxreq = max(maxreq,abs(dqdt(i,k,j,eqn)))
        end do
     end do
  end do

  do i=xs-1,xe+1
  do k=zs-1,ze+1
     q(i,k,0,1)=(8.0d0*uvel(i,k,0)-6.0*q(i,k,1,1)+q(i,k,2,1))*third
     q(i,k,0,2)=(8.0d0*wvel(i,k,0)-6.0*q(i,k,1,2)+q(i,k,2,2))*third
  enddo
  enddo

  write(*,*)'EQN',eqn,'TEST RESIDUAL-Q',maxreq,'PROC',myid

!---------------------------------------------------------------------
  RETURN
END SUBROUTINE RHSQ_test
!*********************************************************************

!*********************************************************************
 SUBROUTINE TEST_PROJECTION

   USE GLOBAL_DATA
   USE mg_solver
   USE MYMPI

   IMPLICIT NONE

   include "parallel.h"

!-----------------------------------------------------------------------
!  Local variables
   INTEGER :: i, k, j, eqn, col(7), error
   INTEGER :: xs,zs,xe,ze,ys,ye,x_end_test
   REAL*8  :: gx,gy,gz,dxa,dya,dza,mdivt2,mdivt1,sum1,sum2
!--------------------------------------------------------------------
!---

   finest_mesh%w => q
   var_updated = 1

   xs = 0
   xe = finest_mesh%xe
   zs = drange(3)
   ze = drange(4)
   ys = drange(5)+1
   ye = drange(6)-1
!
!  the array col contains information about the color of the pints to
! be updated
!  col(1:2), the dimension to be updated, col(3:4), and how many
! point in the
!  j direction need to be updated col(5:6)
!
   col(1)=1; col(2)=2; col(3)=1; col(4) = ndim; 
   col(5)=drange(5); col(6)=drange(6); col(7) = 4

   dxa = 1.0d0/dx
   dza = 1.0d0/dz

   CALL VELOCITY(4)

!
!  evaluate the rhs of the pressure Poisson equation (PPE)
!
   x_end_test = merge(xe+1,xe,is_proc_bnd(2) .AND. ALL(type_bc(1,2,:) == 1))
   sum1=0;sum2=0
   do k=zs,ze
      do i =xs, x_end_test
         do j=ys,ye
            dya = detady(j)/dy
            gx = (uvel(i,k,j)-uvel(i-1,k,j))*dxa
            gy = (vvel(i,k,j)-vvel(i,k,j-1))*dya
            divt(i,k,j)= gx+gy  + rate(i,k,j,4)
!            if(i==x_end_test .AND. is_proc_bnd(2) )print*,i,j,'TEST PROJ',divt(i,k,j),gx,gy,&
!                 'SOURCE',rate(i,k,j,4),'UVEL',&
!                 uvel(i,k,j),uvel(i-1,k,j),vvel(i,k,j),vvel(i,k,j-1),'Q',&
!                 q(i,k,j,1),q(i-1,k,j,1),q(i,k,j,3),q(i,k,j-1,3),&
!                 'RHO',k_p/f(i+1,k,j,1),k_p/f(i,k,j,1),k_p/f(i-1,k,j,1),'MYID is',myid
         enddo
         mdivt1 = maxval(abs(divt(i,k,ys:ye)))
         if( is_proc_bnd(2) ) write(*,*)i,myid,x(i),xloc,'TEST_PROJECTION',mdivt1,dphidx(i,k)
      enddo
   enddo


102 format(i3,1x,a,2f14.6,1x,a,9f14.6)


!-----------------------------------------------------------------------
   RETURN
 END SUBROUTINE TEST_PROJECTION
!*************************************************************************
