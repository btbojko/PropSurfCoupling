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
  REAL*8  :: cplewis,conv,diff,myres(2)
  REAL*8  :: firstcorr,corrflux,nu,s1,s2,s3,sig1,sig2,minmod,g2y
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
  if(dSerialRange(5).eq.0) ys = ys + 1
  if(dSerialRange(6).eq.ny) ye = ye - 1
  zs = drange(3)
  ze = drange(4)
  xs = drange(1)
  xe = drange(2)
  
  myres = zero
  ! Loop over Gas-phase equations
  do eqn=1,neqmax-1
     
     cplewis=cp*lewis(eqn)
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

              nu = abs(vvel(i,k,j))*detady(j)*dy1a
              firstcorr = -gyy/dy2a*nu
              jp = j + merge(0,1,vvel(i,k,j)>0)
              s1 = (f(i,k,jp+1,eqn) - f(i,k,jp,eqn))
              s2 = (f(i,k,jp,eqn) - f(i,k,jp-1,eqn))
              s3 = (f(i,k,jp-1,eqn) - f(i,k,jp-2,eqn))
              if(jp==1) s3=s2;
              if(jp==ny) s1=s2

              sig1 = (sign(half,s1)+sign(half,s2))*min(abs(s1),abs(s2))
              sig2 = (sign(half,s2)+sign(half,s3))*min(abs(s2),abs(s3))

              corrflux = nu * (sig1-sig2)

              gxy = (f(i+1,k,j+1,eqn)-f(i+1,k,j-1,eqn)&
                   -f(i-1,k,j+1,eqn)+f(i-1,k,j-1,eqn))*dxy1
              gzy = (f(i,k+1,j+1,eqn)-f(i,k+1,j-1,eqn)&
                   -f(i,k-1,j+1,eqn)+f(i,k-1,j-1,eqn))*dzy1

              term=1.0d0+dphidx(i,k)**2+dphidz(i,k)**2

              conv = - vvel(i,k,j)*detady(j)*gy &!- firstcorr - corrflux&
                     - uvel(i,k,j)*gx-wvel(i,k,j)*gz
              diff =  lambdag(i,k,j)*(gxx + gzz +&
                   (term)*(&
                   (detady(j)**2)*gyy+deta2dy2(j)*gy)&
                   -2.0d0*dphidx(i,k)*detady(j)*gxy&
                   -2.0d0*dphidz(i,k)*detady(j)*gzy&
                   -(dphi2dx2(i,k)+dphi2dz2(i,k))*detady(j)*gy)/cplewis&
                   +((dlgdx(i,k,j)-dphidx(i,k)*detady(j)*dlgdy(i,k,j))*&
                   (gx-dphidx(i,k)*detady(j)*gy)&
                   +(dlgdz(i,k,j)-dphidz(i,k)*detady(j)*dlgdy(i,k,j))*&
                   (gz-dphidz(i,k)*detady(j)*gy)&
                   +(detady(j)**2)*dlgdy(i,k,j)*gy)/cplewis
               dfdt(i,k,j,eqn) = (conv+diff+rate(i,k,j,eqn))
               myres(1) = max(myres(1),abs(dfdt(i,k,j,eqn)))

! --- value used to evaluate the thermal divergence
!               dfdt(i,k,j,eqn) = conv+diff+rate(i,k,j,eqn)


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

           g2y = (f(i,k,j+1,eqn)-2.0d0*f(i,k,j,eqn)&
                 +f(i,k,j-1,eqn))
           nu = vvel(i,k,0)*detady(j)*dy1a
           firstcorr = -g2y*nu
           jp = j + 1
           s1 = (f(i,k,jp+1,eqn) - f(i,k,jp,eqn))
           s2 = (f(i,k,jp,eqn) - f(i,k,jp-1,eqn))
           s3 = (f(i,k,jp-1,eqn) - f(i,k,jp-2,eqn))
           if(jp==ny) s1=s2

           sig1 = (sign(half,s1)+sign(half,s2))*min(abs(s1),abs(s2))
           sig2 = (sign(half,s2)+sign(half,s3))*min(abs(s2),abs(s3))

           corrflux = nu * (sig1-sig2)

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

           conv = + (vvel(i,k,0)*detady(j)*gy) &!  - firstcorr - corrflux)&
                *rhos(i,k,j)/rhos(i,k,0)
           diff = (gxx + gzz +&
                (term)*(&
                (detady(j)**2)*gyy+deta2dy2(j)*lambdas(i,k,j)*gy)&
                + dphidx(i,k)*detady(j)*gxy&
                + dphidz(i,k)*detady(j)*gzy&
                + (dphi2dx2(i,k)+dphi2dz2(i,k))*detady(j)*&
                lambdas(i,k,j)*gy)/cp
           dfdt(i,k,j,eqn) = conv+diff
           !!dfdt(i,k,j,eqn) = f(i,k,j,eqn) - dfdt(i,k,j,eqn) 
           myres(2) = max(myres(2),abs(dfdt(i,k,j,eqn)))
        end do
     end do
  end do

  call MPI_ALLREDUCE(myres,maxres,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)

!!>  if (myid == 0) write(*,'(i3,a,1x,1p123e16.8)') ncyc,'TEST RESIDUAL',maxres
!  print'(i3,a,1x,1p123e16.8)',myid,'SOLID RES',maxval(abs(dfdt(:,:,:,neqmax)))
  
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
  if(dSerialRange(5).eq.0) ys = ys + 1
  if(dSerialRange(6).eq.ny) ye = ye - 1
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

           conv = - vvel(i,k,j)*gy*ey- firstcorr - corrflux&
                  - uvel(i,k,j)*gx-wvel(i,k,j)*gz
           diff = (dff(i,k,j,1)-dff(i-1,k,j,1))*dx1a&
                 +(dff(i,k,j,2)-dff(i,k-1,j,2))*dz1a&
                 +(dff(i,k,j,3)-dff(i,k,j-1,3))*dyy
!......................................................................
! TO improve efficiency dqdt was previously premultiplied by dtx
!......................................................................
           dqdt(i,k,j,eqn) =  dqdt(i,k,j,eqn)/dtx(i,k,j)+conv+prndtl_cp*diff
           dqdt(i,k,j,eqn) = q(i,k,j,eqn)-dtx(i,k,j)*dqdt(i,k,j,eqn) 
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

!*************************************************
SUBROUTINE CHECK_VEC

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!-----------------------------------------------------------
!!!   Local Variables
  integer i, k, j, n
  REAL*8 :: f_old
  integer ybegin,ytop,ncurb, mynclip
!-----------------------------------------------------------



  ybegin = drange(5)
  ytop = drange(6)

  mynclip = 0
  do j=ybegin,ytop
     do k=drange(3)-1,drange(4)+1
        do i=drange(1)-1,drange(2)+1
           do n=1,neqmax
              if(.NOT. f(i,k,j,n) > min_f(n) .OR. &
                   .NOT. f(i,k,j,n) < max_f(n)) THEN
                 write(*,*)'OUT OF BOUNDS VALUES'
                 write(*,*)i,k,j,n,f(i,k,j,n),min_f(n),max_f(n),&
                      'MYID ID',myid
                 stop 'CHECK_VEC'
              endif
           enddo
        enddo
     enddo
  enddo


!--------------------------------------------------------
  RETURN
END SUBROUTINE CHECK_VEC
!**********************************************************************
!
!*************************************************
SUBROUTINE CHECK_VAL(ff,outcome,neq)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!-----------------------------------------------------------
!!!   Local Variables
  integer i, k, j, n,neq
  REAL*8 :: ff(neq)
  LOGICAL :: outcome
!-----------------------------------------------------------

! TEMP SPECIES BOUNDS MOVED TO MODULE STIFFRATE

!  mx(1) = 4000.0d0
!  mn(1) = 500.0d0
!  do n=2,neqmax-1
!     mn(n) = -1.0
!     mx(n) = 1.d0 - mn(n)
!  enddo


  outcome = .FALSE.
  do n=1,neq
     outcome = outcome .OR. (.NOT. ff(n) > min_f(n)) .OR. (.NOT. ff(n) < max_f(n))
  enddo

!--------------------------------------------------------
  RETURN
END SUBROUTINE CHECK_VAL
!**********************************************************************
!
subroutine CNANV(buffin,L,msg)

  use global_DATA
  implicit none
  
  real*8 :: buffin(L)
  integer :: L
  character*(*) :: msg

  
  if(.not. all(abs(buffin(1:L)) < 1d99)) then
     write(*,*) buffin(1:L)
     write(*,*)myid,'NAN caught'
     write(*,*)msg
     STOP 'CNANV'
  end if

  return
end subroutine CNANV
