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
!     MPI Utilities for ghost value filling 
!
! ********************************************************************

SUBROUTINE FILL_F_SURFACE_GHOST(eqn)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!     Local variables
  INTEGER :: i, k, eqn
  INTEGER :: status(MPI_STATUS_SIZE),xs,zs,xe
  INTEGER ::  ze,xrange,zrange,xsend(2),zsend(2)
!  REAL*8  :: sbuf1(0:maxrange+4),sbuf2(0:maxrange+4),rbuf1(0:maxrange+4),rbuf2(0:maxrange+4)

  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)

  xrange = (xe - xs) + 1
  zrange = (ze - zs) + 1

  xsend(1) = xs
  xsend(2) = xe
  zsend(1) = zs
  zsend(2) = ze


  !!  IF(dSerialRange(1).eq.0)  xsend(1) = 1
  !!  IF(dSerialRange(2).eq.nx) xsend(2) = xe - 1
  !!  IF(dSerialRange(3).eq.0)  zsend(1) = 1
  !!  IF(dSerialRange(4).eq.nz) zsend(2) = ze - 1

  IF(xrange.ne.nx+1) THEN
     DO k = zs,ze
        sbuf1(k) = f(xsend(1),k,0,eqn)
        sbuf2(k) = f(xsend(2),k,0,eqn)
     END DO
     CALL MPI_SENDRECV(sbuf1,zrange,MPI_DOUBLE_PRECISION,&
          dnbr(1),0,rbuf1,zrange,MPI_DOUBLE_PRECISION,&
          dnbr(2),0,comm3d,status,ierr)
     CALL MPI_SENDRECV(sbuf2,zrange,MPI_DOUBLE_PRECISION,&
          dnbr(2),1,rbuf2,zrange,MPI_DOUBLE_PRECISION,&
          dnbr(1),1,comm3d,status,ierr)
     DO k = zs,ze
        f(xs-1,k,0,eqn) = rbuf2(k)
        f(xe+1,k,0,eqn) = rbuf1(k)
     END DO
  ELSE
     DO k=zs,ze
        f(xs-1,k,0,eqn) = f(xe,k,0,eqn)
        f(xe+1,k,0,eqn) = f(xs,k,0,eqn)
     END DO
  END IF

  IF(zrange.ne.nz+1) THEN
     DO k = xs,xe
        sbuf1(k) = f(k,zsend(1),0,eqn)
        sbuf2(k) = f(k,zsend(2),0,eqn)
     END DO
     CALL MPI_SENDRECV(sbuf1,xrange,MPI_DOUBLE_PRECISION,&
          dnbr(3),2,rbuf1,xrange,MPI_DOUBLE_PRECISION,&
          dnbr(4),2,comm3d,status,ierr)
     CALL MPI_SENDRECV(sbuf2,xrange,MPI_DOUBLE_PRECISION,&
          dnbr(4),3,rbuf2,xrange,MPI_DOUBLE_PRECISION,&
          dnbr(3),3,comm3d,status,ierr)
     DO k = xs,xe
        f(k,zs-1,0,eqn) = rbuf2(k)
        f(k,ze+1,0,eqn) = rbuf1(k)
     END DO
  ELSE
     DO i=xs,xe
        f(i,zs-1,0,eqn) = f(i,ze,0,eqn)
        f(i,ze+1,0,eqn) = f(i,zs,0,eqn)
     END DO
  END IF

!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!.................................................................................

  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  IF(issymmetric) THEN
     IF(coords(1) ==  0)           f(-1,:,0,eqn) =   f(1,:,0,eqn)
     IF( coords(1)+1 == pdims(1) ) f(xe+1,:,0,eqn) = f(xe,:,0,eqn)
     IF(dSerialRange(3).eq.0)           f(:,-1,0,eqn) =   f(:,1,0,eqn)
     IF( coords(2)+1 == pdims(2) ) f(:,ze+1,0,eqn) = f(:,ze,0,eqn)
  ENDIF


  RETURN 
END SUBROUTINE FILL_F_SURFACE_GHOST
!*************************************************************************

SUBROUTINE FILL_LAMBDAG_GHOST

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  INTEGER ::  elcount,j,i,zrange,xs,xe,ys,ye,zs,ze,k
!  REAL*8  :: slambda(0:(maxrange+4)**2),rlambda(0:(maxrange+4)**2)
  INTEGER ::  status(MPI_STATUS_SIZE),xrange,yrange


  xrange = (drange(2) - drange(1)) + 1
  yrange = (drange(6) - drange(5)) + 1
  zrange = (drange(4) - drange(3)) + 1

  ys = drange(5)
  ye = drange(6)
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)

!  IF(dSerialRange(1).eq.0) xs = 1
!  IF(dSerialRange(2).eq.nx) xe = xe -1

!  IF(dSerialRange(3).eq.0) zs = 1
!  IF(dSerialRange(4).eq.nz) ze = ze - 1

  IF(xrange.ne.(nx+1)) THEN
     elcount = 0
     DO k = drange(3), drange(4)
        DO j = drange(5), drange(6)
           slambda(elcount) = lambdag(xs,k,j)
           elcount = elcount + 1
        END DO
     END DO
     CALL MPI_SENDRECV(slambda,yrange*zrange,&
          MPI_DOUBLE_PRECISION,dnbr(1),0,rlambda,&
          yrange*zrange,MPI_DOUBLE_PRECISION,dnbr(2),0,&
          comm3d,status,ierr)
     elcount = 0
     DO k = drange(3),drange(4)
        DO j = drange(5), drange(6)
           lambdag(drange(2)+1,k,j) = &
                rlambda(elcount)
           slambda(elcount) = lambdag(xe,k,j)
           elcount = elcount + 1
        END DO
     END DO
     CALL MPI_SENDRECV(slambda,yrange*zrange,&
          MPI_DOUBLE_PRECISION,dnbr(2),1,rlambda,&
          yrange*zrange,MPI_DOUBLE_PRECISION,dnbr(1),1,&
          comm3d,status,ierr)
     elcount = 0
     DO k = drange(3), drange(4)
        DO j = drange(5),drange(6)
           lambdag(drange(1)-1,k,j) = &
                rlambda(elcount)
           elcount = elcount + 1
        END DO
     END DO
  ELSE
     DO k = drange(3),drange(4)
        DO j = ys,ye
           lambdag(xs-1,k,j) = lambdag(xe,k,j)
           lambdag(xe+1,k,j) = lambdag(xs,k,j)
        end do
     END DO
  END IF

  IF(zrange.ne.(nz+1)) THEN
     xs = drange(1) - 1
     xe = drange(2) + 1
     elcount = 0
     DO k = xs, xe
        DO j = ys, ye
           slambda(elcount) = lambdag(k,zs,j)
           elcount = elcount + 1
        END DO
     END DO
     CALL MPI_SENDRECV(slambda,yrange*(xrange+2),&
          MPI_DOUBLE_PRECISION,dnbr(3),2,rlambda,&
          yrange*(xrange+2),MPI_DOUBLE_PRECISION,dnbr(4),2,&
          comm3d,status,ierr)
     elcount = 0
     DO k = xs,xe
        DO j = ys, ye
           lambdag(k,drange(4)+1,j) = rlambda(elcount)
           slambda(elcount) = lambdag(k,ze,j)
           elcount = elcount + 1
        END DO
     END DO

     CALL MPI_SENDRECV(slambda,yrange*(xrange+2),&
          MPI_DOUBLE_PRECISION,dnbr(4),3,rlambda,&
          yrange*(xrange+2),MPI_DOUBLE_PRECISION,dnbr(3),3,&
          comm3d,status,ierr)
     elcount = 0
     DO k = xs,xe
        DO j = ys,ye
           lambdag(k,drange(3)-1,j) = rlambda(elcount)
           elcount = elcount + 1
        END DO
     END DO
  ELSE
     DO i = drange(1)-1,drange(2)+1
        DO j = ys,ye
           lambdag(i,zs-1,j) = lambdag(i,ze,j)
           lambdag(i,ze+1,j) = lambdag(i,zs,j)
        end do
     END DO
  END IF

  IF(yrange.ne.(ny+1)) THEN 
     xs = drange(1)-1
     xe = drange(2)+1
     zs = drange(3)-1
     ze = drange(4)+1
     elcount = 0
     DO i = xs, xe
        DO k = zs, ze
           slambda(elcount) = lambdag(i,k,ye)
           elcount = elcount + 1
        END DO
     END DO
     CALL MPI_SENDRECV(slambda,(xrange+2)*(zrange+2),&
          MPI_DOUBLE_PRECISION,dnbr(6),4,rlambda,&
          (xrange+2)*(zrange+2),MPI_DOUBLE_PRECISION,dnbr(5),4,&
          comm3d,status,ierr) 
     elcount = 0
     DO i = xs, xe
        DO k = zs, ze
           lambdag(i,k,ys-1) = rlambda(elcount)
           slambda(elcount) = lambdag(i,k,ys)
           elcount = elcount + 1
        END DO
     END DO
     CALL MPI_SENDRECV(slambda,(xrange+2)*(zrange+2),&
          MPI_DOUBLE_PRECISION,dnbr(5),5,rlambda,&
          (xrange+2)*(zrange+2),MPI_DOUBLE_PRECISION,dnbr(6),5,&
          comm3d,status,ierr) 
     elcount = elcount + 1
     DO i = xs, xe
        DO k = zs, ze
           lambdag(i,k,ye+1) = rlambda(elcount)
           elcount = elcount + 1
        END DO
     END DO
  END IF

!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!.................................................................................
  IF(issymmetric) THEN
     IF(dSerialRange(1).eq.0)  lambdag(-1,:,:) = lambdag(1,:,:)
     IF(coords(1)+1 == pdims(1) ) lambdag(xe+1,:,:) = lambdag(xe,:,:)
     IF(dSerialRange(3).eq.0)  lambdag(:,-1,:) = lambdag(:,1,:)
     IF(coords(2)+1 == pdims(2) ) lambdag(:,ze+1,:) = lambdag(:,ze,:)
  ENDIF


  RETURN
END SUBROUTINE FILL_LAMBDAG_GHOST

!*************************************************************************
SUBROUTINE FILL_PSI_GHOST

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  INTEGER ::  elcount,j,i,zrange,xs,xe,ys,ye,zs,ze,k
  INTEGER ::  status(MPI_STATUS_SIZE),xrange,yrange,size
!--------------------------------------------------------------  

  xrange = (drange(2) - drange(1)) + 1
  yrange = (drange(6) - drange(5)) + 1
  zrange = (drange(4) - drange(3)) + 1

  ys = drange(5)
  ye = drange(6)
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)

!  IF(dSerialRange(1).eq.0) xs = 1
!  IF(dSerialRange(2).eq.nx) xe = xe -1
!  IF(dSerialRange(3).eq.0) zs = 1
!  IF(dSerialRange(4).eq.nz) ze = ze - 1


  IF(xrange.ne.(nx+1)) THEN
     elcount = 0
     DO k = drange(3), drange(4)
        DO j = drange(5), drange(6)
           slambda(elcount) = psi(xs,k,j)
           elcount = elcount + 1
        END DO
     END DO
     CALL MPI_SENDRECV(slambda,elcount,&
          MPI_DOUBLE_PRECISION,dnbr(1),0,rlambda,&
          elcount,MPI_DOUBLE_PRECISION,dnbr(2),0,&
          comm3d,status,ierr)
     elcount = 0
     DO k = drange(3),drange(4)
        DO j = drange(5), drange(6)
           psi(drange(2)+1,k,j) = &
                rlambda(elcount)
           slambda(elcount) = psi(xe,k,j)
           elcount = elcount + 1
        END DO
     END DO
     CALL MPI_SENDRECV(slambda,yrange*zrange,&
          MPI_DOUBLE_PRECISION,dnbr(2),1,rlambda,&
          yrange*zrange,MPI_DOUBLE_PRECISION,dnbr(1),1,&
          comm3d,status,ierr)
     elcount = 0
     DO k = drange(3), drange(4)
        DO j = drange(5),drange(6)
           psi(drange(1)-1,k,j) = &
                rlambda(elcount)
           elcount = elcount + 1
        END DO
     END DO
  ELSE
     DO k = drange(3),drange(4)
        DO j = ys,ye
           psi(xs-1,k,j) = psi(xe,k,j)
           psi(xe+1,k,j) = psi(xs,k,j)
        end do
     END DO
  END IF

  IF(zrange.ne.(nz+1)) THEN
     xs = drange(1) - 1
     xe = drange(2) + 1
     elcount = 0
     DO k = xs, xe
        DO j = ys, ye
           slambda(elcount) = psi(k,zs,j)
           elcount = elcount + 1
        END DO
     END DO
     CALL MPI_SENDRECV(slambda,yrange*(xrange+2),&
          MPI_DOUBLE_PRECISION,dnbr(3),2,rlambda,&
          yrange*(xrange+2),MPI_DOUBLE_PRECISION,dnbr(4),2,&
          comm3d,status,ierr)
     elcount = 0
     DO k = xs,xe
        DO j = ys, ye
           psi(k,drange(4)+1,j) = rlambda(elcount)
           slambda(elcount) = psi(k,ze,j)
           elcount = elcount + 1
        END DO
     END DO
     CALL MPI_SENDRECV(slambda,yrange*(xrange+2),&
          MPI_DOUBLE_PRECISION,dnbr(4),3,rlambda,&
          yrange*(xrange+2),MPI_DOUBLE_PRECISION,dnbr(3),3,&
          comm3d,status,ierr)
     elcount = 0
     DO k = xs,xe
        DO j = ys,ye
           psi(k,drange(3)-1,j) = rlambda(elcount)
           elcount = elcount + 1
        END DO
     END DO
  ELSE
     DO i = drange(1)-1,drange(2)+1
        DO j = ys,ye
           psi(i,zs-1,j) = psi(i,ze,j)
           psi(i,ze+1,j) = psi(i,zs,j)
        end do
     END DO
  END IF

  IF(yrange.ne.(ny+1)) THEN 
     xs = drange(1)-1
     xe = drange(2)+1
     zs = drange(3)-1
     ze = drange(4)+1
     elcount = 0
     DO i = xs, xe
        DO k = zs, ze
           slambda(elcount) = psi(i,k,ye)
           elcount = elcount + 1
        END DO
     END DO
     CALL MPI_SENDRECV(slambda,(xrange+2)*(zrange+2),&
          MPI_DOUBLE_PRECISION,dnbr(6),4,rlambda,&
          (xrange+2)*(zrange+2),MPI_DOUBLE_PRECISION,dnbr(5),4,&
          comm3d,status,ierr) 
     elcount = 0
     DO i = xs, xe
        DO k = zs, ze
           psi(i,k,ys-1) = rlambda(elcount)
           slambda(elcount) = psi(i,k,ys)
           elcount = elcount + 1
        END DO
     END DO
     CALL MPI_SENDRECV(slambda,(xrange+2)*(zrange+2),&
          MPI_DOUBLE_PRECISION,dnbr(5),5,rlambda,&
          (xrange+2)*(zrange+2),MPI_DOUBLE_PRECISION,dnbr(6),5,&
          comm3d,status,ierr) 
     elcount = elcount + 1
     DO i = xs, xe
        DO k = zs, ze
           psi(i,k,ye+1) = rlambda(elcount)
           elcount = elcount + 1
        END DO
     END DO
  END IF


!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!.................................................................................  
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  IF(issymmetric) THEN
     IF(dSerialRange(1).eq.0)              psi(-1,:,:) =     psi(1,:,:)
     IF(coords(1)+1 == pdims(1) )     psi(xe+1,:,:) =   psi(xe,:,:)
     IF(dSerialRange(3).eq.0)              psi(:,-1,:) =     psi(:,1,:)
     IF(coords(2)+1 == pdims(2) )     psi(:,ze+1,:) =   psi(:,ze,:)
  ENDIF

!-----------------------------------------------------
  RETURN
END SUBROUTINE FILL_PSI_GHOST
!********************************************************************

SUBROUTINE FILL_RB_GHOST

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!     Local Variables
  INTEGER ::  i, k,elcount,xs,xe,zs,ze
  INTEGER ::  xrange,zrange,status(MPI_STATUS_SIZE)
  REAL*8  :: rbs_left(0:(3*(maxrange+4))),rbr_left(0:(3*(maxrange+4)))
  REAL*8  :: rbs_right(0:(3*(maxrange+4))),rbr_right(0:(3*(maxrange+4)))

  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)


  xrange = (xe - xs)+1
  zrange = (ze - zs)+1

  IF(xrange.ne.(nx+1)) THEN
!     IF(dSerialRange(1).eq.0) xs = 1
!     IF(dSerialRange(2).eq.nx) xe = xe - 1
     elcount = 0
     DO i = zs,ze
        rbs_left(elcount) = rb(xs,i)
        rbs_right(elcount) = rb(xe,i)
        elcount = elcount + 1
     END DO
     CALL MPI_SENDRECV(rbs_left,zrange,&
          MPI_DOUBLE_PRECISION,dnbr(1),0,rbr_right,&
          zrange,MPI_DOUBLE_PRECISION,dnbr(2),0,&
          comm3d,status,ierr)
     CALL MPI_SENDRECV(rbs_right,zrange,&
          MPI_DOUBLE_PRECISION,dnbr(2),1,rbr_left,&
          zrange,MPI_DOUBLE_PRECISION,dnbr(1),1,&
          comm3d,status,ierr)
     elcount = 0
     DO i = zs,ze
        rb(drange(2)+1,i) = rbr_right(elcount)
        rb(drange(1)-1,i) = rbr_left(elcount)
        rbs_left(elcount) = rb(xs+1,i)
        rbs_right(elcount) = rb(xe-1,i)
        elcount = elcount + 1
     END DO
     CALL MPI_SENDRECV(rbs_left,zrange,&
          MPI_DOUBLE_PRECISION,dnbr(1),2,rbr_right,&
          zrange,MPI_DOUBLE_PRECISION,dnbr(2),2,&
          comm3d,status,ierr)
     CALL MPI_SENDRECV(rbs_right,zrange,&
          MPI_DOUBLE_PRECISION,dnbr(2),3,rbr_left,&
          zrange,MPI_DOUBLE_PRECISION,dnbr(1),3,&
          comm3d,status,ierr)
     elcount = 0
     DO i = zs,ze
        rb(drange(2)+2,i) = rbr_right(elcount)
        rb(drange(1)-2,i) = rbr_left(elcount)
        elcount = elcount + 1
     END DO
  ELSE
     DO k=drange(3),drange(4)
        rb(xs-1,k) = rb(xe,k)
        rb(xs-2,k) = rb(xe-1,k)
        rb(xe+1,k) = rb(xs,k)
        rb(xe+2,k) = rb(xs+1,k)
     END DO
  END IF
  xs = drange(1)
  xe = drange(2)
  IF(zrange.ne.(nz+1)) THEN
!     IF(dSerialRange(3).eq.0) zs = 1
!     IF(dSerialRange(4).eq.nz) ze = ze - 1
     elcount = 0
     DO i = xs,xe
        rbs_left(elcount) = rb(i,zs)
        rbs_right(elcount) = rb(i,ze)
        elcount = elcount + 1
     END DO
     CALL MPI_SENDRECV(rbs_left,xrange,&
          MPI_DOUBLE_PRECISION,dnbr(3),4,rbr_right,&
          xrange,MPI_DOUBLE_PRECISION,dnbr(4),4,&
          comm3d,status,ierr)
     CALL MPI_SENDRECV(rbs_right,xrange,&
          MPI_DOUBLE_PRECISION,dnbr(4),5,rbr_left,&
          xrange,MPI_DOUBLE_PRECISION,dnbr(3),5,&
          comm3d,status,ierr)
     elcount = 0
     DO i = xs,xe
        rb(i,drange(3)-1) = rbr_left(elcount)
        rb(i,drange(4)+1) = rbr_right(elcount)
        rbs_left(elcount) = rb(i,zs+1)
        rbs_right(elcount) = rb(i,ze-1)
        elcount = elcount + 1
     END DO
     CALL MPI_SENDRECV(rbs_left,xrange,&
          MPI_DOUBLE_PRECISION,dnbr(3),6,rbr_right,&
          xrange+4,MPI_DOUBLE_PRECISION,dnbr(4),6,&
          comm3d,status,ierr)
     CALL MPI_SENDRECV(rbs_right,xrange,&
          MPI_DOUBLE_PRECISION,dnbr(4),7,rbr_left,&
          xrange,MPI_DOUBLE_PRECISION,dnbr(3),7,&
          comm3d,status,ierr)
     elcount = 0
     DO i = xs,xe
        rb(i,drange(3)-2) = rbr_left(elcount)
        rb(i,drange(4)+2) = rbr_right(elcount)
        elcount = elcount + 1
     END DO
  ELSE
     DO i=drange(1),drange(2)
        rb(i,zs-1) = rb(i,ze)
        rb(i,zs-2) = rb(i,ze-1)
        rb(i,ze+1) = rb(i,zs)
        rb(i,ze+2) = rb(i,zs+1)
     END DO
  END IF


!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!.................................................................................  
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  IF(issymmetric) THEN
     IF(dSerialRange(1).eq.0)       rb(-1,:) =      rb(1,:)
     IF(dSerialRange(1).eq.0)       rb(-2,:) =      rb(2,:)
     IF(coords(1)+1 == pdims(1) ) THEN
        rb(xe+1,:) =    rb(xe,:)
        rb(xe+2,:) =    rb(xe-1,:)
     ENDIF
     IF(dSerialRange(3).eq.0)       rb(:,-1) =      rb(:,1)
     IF(dSerialRange(3).eq.0)       rb(:,-2) =      rb(:,2)
     IF(coords(2)+1 == pdims(2) ) THEN
        rb(:,ze+1) =    rb(:,ze)
        rb(:,ze+2) =    rb(:,ze-1)
     ENDIF
  ENDIF

!----------------------------------------------------------
  RETURN
END SUBROUTINE FILL_RB_GHOST

!*************************************************************
SUBROUTINE FILL_F_GHOST(eqn)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!     Local variables
  INTEGER ::  i, j, k, eqn
  INTEGER ::  status(MPI_STATUS_SIZE),yrange,xrange,zrange
  INTEGER ::  ys,ye,xs,xe,zs,ze,elcount
!  REAL*8  :: myright_ff_sbuffer(0:(xr(3)+2)**2),myright_ff_rbuffer(0:(xr(3)+2)**2)
!  REAL*8  :: myleft_ff_sbuffer(0:(xr(3)+2)**2),myleft_ff_rbuffer(0:(xr(3)+2)**2)

  ys = drange(5)
  ye = drange(6)
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)

  yrange = (ye-ys)+1
  xrange = (xe-xs)+1
  zrange = (ze-zs)+1


  IF(xrange.ne.(nx+1)) THEN
!     IF(dSerialRange(1).eq.0) xs = 1
!     IF(dSerialRange(2).eq.nx) xe = xe - 1
     elcount = 0
     DO k = zs, ze
        DO j = ys, ye
           myleft_ff_sbuffer(elcount)=f(xs,k,j,eqn)
           myright_ff_sbuffer(elcount)=f(xe,k,j,eqn)
           elcount = elcount + 1
        END DO
     END DO
!
     CALL MPI_SENDRECV(myleft_ff_sbuffer(0),yrange*zrange,&
          MPI_DOUBLE_PRECISION,&
          dnbr(1),0,myright_ff_rbuffer(0),yrange*zrange,&
          MPI_DOUBLE_PRECISION,&
          dnbr(2),0,comm3d,status,ierr)
     CALL MPI_SENDRECV(myright_ff_sbuffer(0),yrange*zrange,&
          MPI_DOUBLE_PRECISION,&
          dnbr(2),1,myleft_ff_rbuffer(0),yrange*zrange,&
          MPI_DOUBLE_PRECISION,&
          dnbr(1),1,comm3d,status,ierr)
!
     xs = drange(1)
     xe = drange(2)
     elcount = 0
     DO k = zs,ze
        DO j = ys, ye
           f(xe+1,k,j,eqn)=myright_ff_rbuffer(elcount)
           f(xs-1,k,j,eqn)=myleft_ff_rbuffer(elcount)
           elcount = elcount + 1
        END DO
     END DO
  ELSE
     DO k = zs,ze
        DO j = ys,ye
           f(xe+1,k,j,eqn) = f(xs,k,j,eqn)
           f(xs-1,k,j,eqn) = f(xe,k,j,eqn)
        END DO
     END DO
  END IF

  IF(zrange.ne.(nz+1)) THEN
!     IF(dSerialRange(3).eq.0) zs = 1
!     IF(dSerialRange(4).eq.nz) ze = ze - 1
     xs = drange(1) 
     xe = drange(2) 
     elcount = 0
     DO i = xs-1, xe+1
        DO j = ys, ye
           myleft_ff_sbuffer(elcount)=f(i,zs,j,eqn)
           myright_ff_sbuffer(elcount)=f(i,ze,j,eqn)
           elcount = elcount + 1
        END DO
     END DO
     CALL MPI_SENDRECV(myleft_ff_sbuffer,elcount,&
          MPI_DOUBLE_PRECISION,&
          dnbr(3),2,myright_ff_rbuffer,elcount,&
          MPI_DOUBLE_PRECISION,&
          dnbr(4),2,comm3d,status,ierr)
     CALL MPI_SENDRECV(myright_ff_sbuffer,elcount,&
          MPI_DOUBLE_PRECISION,&
          dnbr(4),3,myleft_ff_rbuffer,elcount,&
          MPI_DOUBLE_PRECISION,&
          dnbr(3),3,comm3d,status,ierr)
     zs = drange(3)
     ze = drange(4)
     elcount = 0
     DO i = xs-1,xe+1
        DO j = ys, ye
           f(i,ze+1,j,eqn)=myright_ff_rbuffer(elcount)
           f(i,zs-1,j,eqn)=myleft_ff_rbuffer(elcount)
           elcount = elcount + 1
        END DO
     END DO
     xs = drange(1)
     xe = drange(2)
  ELSE
     xs = drange(1) 
     xe = drange(2)
     DO i = xs-1,xe+1
        DO j = ys,ye
           f(i,ze+1,j,eqn) = f(i,zs,j,eqn)
           f(i,zs-1,j,eqn) = f(i,ze,j,eqn)
        END DO
     END DO
     xs = drange(1)
     xe = drange(2)
  END IF

  IF(yrange.ne.(ny+1)) THEN
     write(*,*)'DOING a swap in y ??'
  END IF
!verbcorn

!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!.................................................................................  
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  IF(issymmetric) THEN
     IF(dSerialRange(1).eq.0)          f(-1,:,:,eqn) =   f(1,:,:,eqn)
     IF(coords(1)+1 == pdims(1) ) f(xe+1,:,:,eqn) = f(xe,:,:,eqn)
     IF(dSerialRange(3).eq.0)          f(:,-1,:,eqn) =   f(:,1,:,eqn)
     IF(coords(2)+1 == pdims(2) ) f(:,ze+1,:,eqn) = f(:,ze,:,eqn)
  ENDIF


  RETURN
END SUBROUTINE FILL_F_GHOST

!
!
! ********************************************************************
!
!
SUBROUTINE FILL_PHI_GHOST
  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

  INTEGER ::  i,k,xrange,zrange,yrange,xs,xe,zs,ze,ys,ye
  INTEGER ::  status(MPI_STATUS_SIZE),xsend(2),zsend(2),elcount
  REAL*8  :: sbuf_left(0:maxrange+4),sbuf_right(0:maxrange+4),sbuf_left2(0:maxrange+4),sbuf_right2(0:maxrange+4)
  REAL*8  :: rbuf_right(0:maxrange+4),rbuf_left2(0:maxrange+4),rbuf_left(0:maxrange+4),rbuf_right2(0:maxrange+4)

  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  ys = drange(5)
  ye = drange(6)

  xrange = (xe - xs) + 1
  zrange = (ze - zs) + 1
  yrange = (ye - ys) + 1

  IF(xrange.ne.nx+1) THEN
     xsend(1) = xs
     xsend(2) = xe
!     IF(dSerialRange(1).eq.0) xsend(1) = 1
!     IF(dSerialRange(2).eq.nx) xsend(2) = xe -1
     elcount = 0
     DO k = zs,ze
        sbuf_left(elcount) = phi(xsend(1),k)
        sbuf_right(elcount) = phi(xsend(2),k)
        sbuf_left2(elcount) = phi(xsend(1)+1,k)
        sbuf_right2(elcount) = phi(xsend(2)-1,k)
        elcount = elcount + 1
     END DO

     CALL MPI_SENDRECV(sbuf_left,zrange,MPI_DOUBLE_PRECISION,&
          dnbr(1),0,rbuf_right,zrange,MPI_DOUBLE_PRECISION,dnbr(2),&
          0,comm3d,status,ierr)
     CALL MPI_SENDRECV(sbuf_left2,zrange,MPI_DOUBLE_PRECISION,&
          dnbr(1),1,rbuf_right2,zrange,MPI_DOUBLE_PRECISION,dnbr(2),&
          1,comm3d,status,ierr)
     CALL MPI_SENDRECV(sbuf_right,zrange,MPI_DOUBLE_PRECISION,&
          dnbr(2),2,rbuf_left,zrange,MPI_DOUBLE_PRECISION,dnbr(1),&
          2,comm3d,status,ierr)
     CALL MPI_SENDRECV(sbuf_right2,zrange,MPI_DOUBLE_PRECISION,&
          dnbr(2),3,rbuf_left2,zrange,MPI_DOUBLE_PRECISION,dnbr(1),&
          3,comm3d,status,ierr)
     elcount = 0
     DO k = zs,ze
        phi(drange(1)-1,k) = rbuf_left(elcount)
        phi(drange(1)-2,k) = rbuf_left2(elcount)
        phi(drange(2)+1,k) = rbuf_right(elcount)
        phi(drange(2)+2,k) = rbuf_right2(elcount)
        elcount = elcount + 1
     END DO
  ELSE
     DO k=zs,ze
        phi(xs-1,k) = phi(xe,k)
        phi(xs-2,k) = phi(xe-1,k)
        phi(xe+1,k) = phi(xs,k)
        phi(xe+2,k) = phi(xs+1,k)
     END DO
  END IF

  IF(zrange.ne.nz+1) THEN
     xs = drange(1)
     xe = drange(2)
     zs = drange(3)
     ze = drange(4)
     ys = drange(5)
     ye = drange(6)
     xrange = (xe - xs) + 1
     zrange = (ze - zs) + 1
     yrange = (ye - ys) + 1
     zsend(1) = zs
     zsend(2) = ze
!     IF(dSerialRange(3).eq.0) zsend(1) = 1
!     IF(dSerialRange(4).eq.nz) zsend(2) = ze - 1
     elcount = 0
     DO i = xs-1, xe+1
        sbuf_left(elcount) = phi(i,zsend(1))
        sbuf_left2(elcount) = phi(i,zsend(1)+1)
        sbuf_right(elcount) = phi(i,zsend(2))
        sbuf_right2(elcount) = phi(i,zsend(2)-1)
        elcount = elcount + 1
     END DO
!
     CALL MPI_SENDRECV(sbuf_left,elcount,MPI_DOUBLE_PRECISION,&
          dnbr(3),4,rbuf_right,elcount,MPI_DOUBLE_PRECISION,&
          dnbr(4),4,comm3d,status,ierr)
     CALL MPI_SENDRECV(sbuf_left2,elcount,MPI_DOUBLE_PRECISION,&
          dnbr(3),5,rbuf_right2,elcount,MPI_DOUBLE_PRECISION,&
          dnbr(4),5,comm3d,status,ierr)
     CALL MPI_SENDRECV(sbuf_right,elcount,MPI_DOUBLE_PRECISION,&
          dnbr(4),6,rbuf_left,elcount,MPI_DOUBLE_PRECISION,&
          dnbr(3),6,comm3d,status,ierr)
     CALL MPI_SENDRECV(sbuf_right2,elcount,MPI_DOUBLE_PRECISION,&
          dnbr(4),7,rbuf_left2,elcount,MPI_DOUBLE_PRECISION,&
          dnbr(3),7,comm3d,status,ierr)
     elcount = 0
     DO i = xs-1, xe+1
        phi(i,ze+1) = rbuf_right(elcount)
        phi(i,ze+2) = rbuf_right2(elcount)
        phi(i,zs-1) = rbuf_left(elcount)
        phi(i,zs-2) = rbuf_left2(elcount)
        elcount = elcount + 1
     END DO
!
  ELSE
     DO i=xs-1,xe+1
        phi(i,zs-1) = phi(i,ze)
        phi(i,zs-2) = phi(i,ze-1)
        phi(i,ze+1) = phi(i,zs)
        phi(i,ze+2) = phi(i,zs+1)
     END DO
  END IF

!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!.................................................................................  
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  IF(issymmetric) THEN
     IF(dSerialRange(1).eq.0)      phi(-1,:) =     phi(1,:)
     IF(dSerialRange(1).eq.0)      phi(-2,:) =     phi(2,:)
     IF(coords(1)+1 == pdims(1) ) THEN
        phi(xe+1,:) =   phi(xe,:)
        phi(xe+2,:) =   phi(xe-1,:)
     ENDIF
     IF(dSerialRange(3).eq.0)      phi(:,-1) =     phi(:,1)
     IF(dSerialRange(3).eq.0)      phi(:,-2) =     phi(:,2)
     IF(coords(2)+1 == pdims(2) ) THEN
        phi(:,ze+1) =   phi(:,ze)
        phi(:,ze+2) =   phi(:,ze-1)
     ENDIF
  ENDIF

  RETURN
END SUBROUTINE FILL_PHI_GHOST

!---------------------------------------------------
SUBROUTINE FILL_FNEW

  USE GLOBAL_DATA
  USE BC_DATA
  USE MYMPI
  IMPLICIT NONE    

!     Local variables
  INTEGER ::  i, k,status(MPI_STATUS_SIZE),xs,zs,xe
  INTEGER ::  ze,xrange,zrange,xsend(2),zsend(2)
!  REAL*8  :: sbuf1(0:maxrange+4),sbuf2(0:maxrange+4),rbuf1(0:maxrange+4),rbuf2(0:maxrange+4)

  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)

  xrange = (xe - xs) + 1
  zrange = (ze - zs) + 1

  xsend(1) = xs
  xsend(2) = xe
  zsend(1) = zs
  zsend(2) = ze


!  IF(dSerialRange(1).eq.0)  xsend(1) = 1
!  IF(dSerialRange(2).eq.nx) xsend(2) = xe - 1
!  IF(dSerialRange(3).eq.0)  zsend(1) = 1
!  IF(dSerialRange(4).eq.nz) zsend(2) = ze - 1


  IF(xrange.ne.nx+1) THEN
     DO k = zs,ze
        sbuf1(k) = ffnew(xsend(1),k)
        sbuf2(k) = ffnew(xsend(2),k)
     END DO
     CALL MPI_SENDRECV(sbuf1,zrange,MPI_DOUBLE_PRECISION,&
          dnbr(1),0,rbuf1,zrange,MPI_DOUBLE_PRECISION,&
          dnbr(2),0,comm3d,status,ierr)
     CALL MPI_SENDRECV(sbuf2,zrange,MPI_DOUBLE_PRECISION,&
          dnbr(2),1,rbuf2,zrange,MPI_DOUBLE_PRECISION,&
          dnbr(1),1,comm3d,status,ierr)
     DO k = zs,ze
        ffnew(xs-1,k) = rbuf2(k)
        ffnew(xe+1,k) = rbuf1(k)
     END DO
  ELSE
     DO k=zs,ze
        ffnew(xs-1,k) = ffnew(xe,k)
        ffnew(xe+1,k) = ffnew(xs,k)
     END DO
  END IF
  IF(zrange.ne.nz+1) THEN
     DO k = xs,xe
        sbuf1(k) = ffnew(k,zsend(1))
        sbuf2(k) = ffnew(k,zsend(2))
     END DO
     CALL MPI_SENDRECV(sbuf1,xrange,MPI_DOUBLE_PRECISION,&
          dnbr(3),2,rbuf1,xrange,MPI_DOUBLE_PRECISION,&
          dnbr(4),2,comm3d,status,ierr)
     CALL MPI_SENDRECV(sbuf2,xrange,MPI_DOUBLE_PRECISION,&
          dnbr(4),3,rbuf2,xrange,MPI_DOUBLE_PRECISION,&
          dnbr(3),3,comm3d,status,ierr)
     DO k = xs,xe
        ffnew(k,zs-1) = rbuf2(k)
        ffnew(k,ze+1) = rbuf1(k)
     END DO
  ELSE
     DO i=xs,xe
        ffnew(i,zs-1) = ffnew(i,ze)
        ffnew(i,ze+1) = ffnew(i,zs)
     END DO
  END IF
!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!.................................................................................  
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  IF(issymmetric) THEN
     IF(dSerialRange(1).eq.0)          ffnew(-1,:) =   ffnew(1,:)
     IF(coords(1)+1 == pdims(1) ) ffnew(xe+1,:) = ffnew(xe,:)
     IF(dSerialRange(3).eq.0)          ffnew(:,-1) =   ffnew(:,1)
     IF(coords(2)+1 == pdims(2) ) ffnew(:,ze+1) = ffnew(:,ze)
  ENDIF


!------------------------------------------------------------------------------------
  RETURN
END SUBROUTINE FILL_FNEW
!*************************************************************************************



!*************************************************************************************
SUBROUTINE parallel_swap(px,col,current_mesh)

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE    

!------------------------------------------------------------------------------------------------
! Dummy variables:
  REAL(KIND=double), DIMENSION(:,:,:), POINTER  :: px
  TYPE(mesh_pointers),INTENT(IN) :: current_mesh

  INTEGER, INTENT(IN) :: col(6)

! Local variables
  INTEGER ::  elcount,j,i,zrange,xs,xe,ys,ye,zs,ze,k,iadd,kadd
  INTEGER ::  status(MPI_STATUS_SIZE),tmpdnbr(4),tmpcomm
!------------------------------------------------------------------------------------------------

  xs = 0
  xe = col(2)
  zs = 0
  ze = col(4)
  iadd=0;kadd=0;
  tmpdnbr(1:4) = current_mesh%dnbr(1:4)
  tmpcomm = current_mesh%comm3d




  IF(isparallel(1)) THEN      !parallel in x

     elcount = 0
     DO k = col(3), col(4)
        DO j = col(5), col(6)
           slambda(elcount) = px(xe,k+kadd,j)    !SEND the last colon, kadd
           elcount = elcount + 1
        END DO
     END DO

     CALL MPI_SENDRECV(slambda,elcount,&
          MPI_DOUBLE_PRECISION,tmpdnbr(2),0,rlambda,&
          elcount,MPI_DOUBLE_PRECISION,tmpdnbr(1),0,&
          tmpcomm,status,ierr)

     elcount = 0
     DO k = col(3), col(4)
        DO j = col(5), col(6)
           px(xs-1,k+kadd,j) = rlambda(elcount)   !RECEIVE the last colon, kadd
           slambda(elcount) = px(xs,k,j)          !SEND the first colon
           elcount = elcount + 1
        END DO
     END DO

     CALL MPI_SENDRECV(slambda,elcount,&
          MPI_DOUBLE_PRECISION,tmpdnbr(1),1,rlambda,&
          elcount,MPI_DOUBLE_PRECISION,tmpdnbr(2),1,&
          tmpcomm,status,ierr)

     elcount = 0
     DO k = col(3), col(4)
        DO j = col(5), col(6)
           px(xe+1,k,j) = rlambda(elcount)    !RECEIVE the first colon
           elcount = elcount + 1
        END DO
     END DO


  ELSE
     DO k = col(3), col(4)
        DO j = col(5), col(6)
           px(xs-1,k+kadd,j) = px(xe,k+kadd,j)
           px(xe+1,k,j) = px(xs,k,j)
        END DO
     END DO
  END IF    !isparallel(1) x direction

  IF(isparallel(2)) THEN    !parallel in z direction
     elcount = 0
     DO i = col(1)-1, col(2)+1
        DO j = col(5), col(6)
           slambda(elcount) = px(i+iadd,ze,j)
           elcount = elcount + 1
        END DO
     END DO
     CALL MPI_SENDRECV(slambda,elcount,&
          MPI_DOUBLE_PRECISION,tmpdnbr(4),2,rlambda,&
          elcount,MPI_DOUBLE_PRECISION,tmpdnbr(3),2,&
          tmpcomm,status,ierr)
     elcount = 0
     DO i = col(1)-1, col(2)+1
        DO j = col(5), col(6)
           px (i+iadd,zs-1,j) = rlambda(elcount)
           slambda(elcount) = px(i,zs,j)
           elcount = elcount + 1
        END DO
     END DO

     CALL MPI_SENDRECV(slambda,elcount,&
          MPI_DOUBLE_PRECISION,tmpdnbr(3),3,rlambda,&
          elcount,MPI_DOUBLE_PRECISION,tmpdnbr(4),3,&
          tmpcomm,status,ierr)

     elcount = 0
     DO i = col(1)-1, col(2)+1
        DO j = col(5), col(6)
           px(i,ze+1,j) = rlambda(elcount)
           elcount = elcount + 1
        END DO
     END DO

  ELSE

     DO i = col(1)-1, col(2)+1
        DO j = col(5), col(6)
           px (i+iadd,zs-1,j) = px (i+iadd,ze,j)
           px (i,ze+1,j) = px (i,zs,j)
        end do
     END DO

  END IF

  IF(isparallel(3)) THEN 
     write(*,*) myid,':: Doing a Y send?!?'
  END IF


!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!.................................................................................
  IF(issymmetric) THEN
     IF(dSerialRange(1).eq.0)         px(-1,:,:) =   px(1,:,:)
     IF(coords(1)+1 == pdims(1)) px(xe+1,:,:) = px(xe,:,:)
     IF(dSerialRange(3).eq.0)         px(:,-1,:) =   px(:,1,:)
     IF(coords(2)+1 == pdims(2)) px(:,ze+1,:) = px(:,ze,:)
  ENDIF



!------------------------------------------------------------------------
  RETURN
END SUBROUTINE parallel_swap
!*************************************************************************


!*************************************************************************
SUBROUTINE parallel_swap4(col,current_mesh)

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE    

!------------------------------------------------------------------------
! Dummy variables:
  TYPE(mesh_pointers),INTENT(IN) :: current_mesh
  INTEGER, INTENT(IN) :: col(7)

! Local variables
  REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: pw

  INTEGER ::  i,k,j,xs,xe,ys,ye,zs,ze,icolor,color_start,color_end, indx, iss(maxsize)
  INTEGER ::  status(MPI_STATUS_SIZE),xrange,yrange,size,eqn,eqn1,neqn,maxmem
  INTEGER ::  indsend(4),indrecv(4), nsend,nrecv, snbr,rnbr, is,ks, elcount, elrecv
  INTEGER ::  iindrecv(4),ii,kk,i1,i2,k1,k2
  INTEGER ::  tmpdnbr(4),tmpcomm
  REAL(KIND=double) ::  ss(maxsize)
!--------------------------------------------------------------------------

  if (col(1) > 0) then
     pw =>  current_mesh%w; 
  else
     pw => wrk4
  endif
!
  color_start= abs(col(1))
  color_end= col(2)
  eqn1 = col (3)
  neqn = col (4)
  ys = col(5)
  ye = col(6)
  xe = current_mesh%xe ; ze = current_mesh%ze;
  tmpdnbr(1:4) = current_mesh%dnbr(1:4)
  tmpcomm = current_mesh%comm3d
  indsend(1)=0;indsend(2)=xe;indsend(3)=0;indsend(4)=ze;
  indrecv(1)=-1;indrecv(2)=xe+1;indrecv(3)=-1;indrecv(4)=ze+1;

  iindrecv(1)=-2;iindrecv(2)=xe+2;iindrecv(3)=-2;iindrecv(4)=ze+2;

10 continue

  maxmem=min(minval(shape(slambda)),minval(shape(rlambda)))


  IF(isparallel(1)) THEN      !parallel in x

     neighloop_x: DO snbr =2,1,-1
        rnbr=rmap(snbr)

        i=indsend(snbr)
        elcount = 0
        do icolor = color_start,color_end 
           nsend= current_mesh%nsnd(snbr,icolor)
           do indx = 1,nsend
              k= current_mesh%sndclr(indx,snbr,icolor)
              do j = ys,ye
                 do eqn = eqn1,neqn
                    slambda(elcount) = pw(i,k,j,eqn)
                    elcount=elcount+1
                 enddo
              enddo
           enddo
        enddo
!
! check the dimension of the buffer used in the communication
!
        if(elcount.gt.maxmem) THEN 
           write(*,*) 'INSUFFICIENT MEMORY IN COMMUNICATION',elcount,maxmem
           DEALLOCATE(slambda,rlambda)
           ALLOCATE(slambda(0:elcount+1),rlambda(0:elcount+1))
           GOTO 10
        endif
        elrecv=sum(current_mesh%nrcv(rnbr,color_start:color_end))*(ye-ys+1)*(neqn-eqn1+1)
        CALL MPI_SENDRECV(slambda,elcount,&
             MPI_DOUBLE_PRECISION,dnbr(snbr),0,rlambda,&
             elrecv,MPI_DOUBLE_PRECISION,dnbr(rnbr),0,&
             tmpcomm,status,ierr)

        i=indrecv(rnbr)
        ii=iindrecv(rnbr)
        i1 = min(i,ii); i2 = max(i,ii)
        elcount=0
        do icolor = color_start,color_end 
           nrecv= current_mesh%nrcv(rnbr,icolor)
           do indx = 1,nrecv
              k= current_mesh%rcvclr(indx,rnbr,icolor)
              do j = ys,ye
                 do eqn = eqn1,neqn
                    pw(i1:i2,k,j,eqn) = rlambda(elcount)
                    elcount=elcount+1
                 enddo
              enddo
           enddo
        enddo

     ENDDO neighloop_x

  ELSE
!
     DO snbr =1,2
        rnbr=rmap(snbr)
        i=indrecv(rnbr);is=indsend(snbr);ii=iindrecv(rnbr);
        i1 = min(i,ii); i2 = max(i,ii)
        do icolor = color_start,color_end
           nrecv= current_mesh%nrcv(rnbr,icolor)
           do indx = 1,nrecv
              k= current_mesh%rcvclr(indx,rnbr,icolor)
              do j = ys,ye
                 do eqn = eqn1,neqn
                    pw(i1:i2,k,j,eqn)=pw(is,k,j,eqn)
                 enddo
              enddo
           enddo
        enddo
     ENDDO

  ENDIF

  IF(isparallel(2)) THEN   !parallel in z

     DO snbr = 3,4
        rnbr=rmap(snbr)

        k=indsend(snbr)
        elcount = 0
        do icolor = color_start,color_end 
           nsend= current_mesh%nsnd(snbr,icolor)
           do indx = 1,nsend
              i= current_mesh%sndclr(indx,snbr,icolor)
              do j = ys,ye
                 do eqn = eqn1,neqn
                    slambda(elcount) = pw(i,k,j,eqn)
                    elcount=elcount+1
                 enddo
              enddo
           enddo
        enddo

!
! check the dimension of the buffer used in the communication
!      
        if(elcount.gt.maxmem) THEN 
           write(*,*) 'INSUFFICIENT MEMORY IN COMMUNICATION',elcount,maxmem
           DEALLOCATE(slambda,rlambda)
           ALLOCATE(slambda(0:elcount+1),rlambda(0:elcount+1))
           GOTO 10
        endif

        elrecv=sum(current_mesh%nrcv(rnbr,color_start:color_end))*(ye-ys+1)*(neqn-eqn1+1)
        CALL MPI_SENDRECV(slambda,elcount,&
             MPI_DOUBLE_PRECISION,dnbr(snbr),0,rlambda,&
             elrecv,MPI_DOUBLE_PRECISION,dnbr(rnbr),0,&
             tmpcomm,status,ierr)

        k=indrecv(rnbr);kk = iindrecv(rnbr)
        k1 = min(k,kk);k2 = max(k,kk)
        elcount=0
        do icolor = color_start,color_end 
           nrecv= current_mesh%nrcv(rnbr,icolor)
           do indx = 1,nrecv
              i= current_mesh%rcvclr(indx,rnbr,icolor)
              do j = ys,ye
                 do eqn = eqn1,neqn
                    pw(i,k1:k2,j,eqn) = rlambda(elcount)
                    elcount=elcount+1
                 enddo
              enddo
           enddo
        enddo

     ENDDO

  ELSE
!
     DO snbr =3,4
        rnbr=rmap(snbr)
        k=indrecv(rnbr);ks=indsend(snbr);kk = iindrecv(rnbr)
        k1 = min(k,kk);k2 = max(k,kk)
        do icolor = color_start,color_end
           nrecv= current_mesh%nrcv(rnbr,icolor)
           do indx = 1,nrecv
              i= current_mesh%rcvclr(indx,rnbr,icolor)
              do j = ys,ye
                 do eqn = eqn1,neqn
                    pw(i,k1:k2,j,eqn)=pw(i,ks,j,eqn)
                 enddo
              enddo
           enddo
        enddo
     ENDDO

  ENDIF


!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!.................................................................................
  IF(issymmetric) THEN
     DO eqn = eqn1,neqn

!     ss=min(dble(2*eqn) - 3.0d0 , one )    !eqn=1,2,3  ss=-1,1,1   iss=-1,0,0
!     iss = min(2*eqn - 3 , 0)

        ss = one
        iss = 0
        IF(col(7) /= 1) THEN 
           ss(1) = -1.0d0
           iss(1) = -1
        ENDIF

        IF(dSerialRange(1).eq.0)           pw(-1,:,:,eqn) =   ss(eqn) * pw(1+iss(eqn),:,:,eqn)
        IF(coords(1)+1 == pdims(1) )  pw(xe+1,:,:,eqn) = ss(eqn) * pw(xe+iss(eqn),:,:,eqn)


        IF(dSerialRange(3).eq.0)          pw(:,-1,:,eqn) =   ss(eqn) * pw(:,1+iss(eqn),:,eqn)
        IF(coords(2)+1 == pdims(2))  pw(:,ze+1,:,eqn) = ss(eqn) * pw(:,ze+iss(eqn),:,eqn)

!---2D only
        if(ipack == 0) then
           pw(:,-1,:,eqn) =   pw(:,0,:,eqn)
           pw(:,ze+1,:,eqn) = pw(:,ze,:,eqn)
        end if
     ENDDO
  ENDIF

!------------------------------------------------------------------------
  RETURN
END SUBROUTINE parallel_swap4
!************************************************************************

!**********************************************************************************
SUBROUTINE color_swap(current_mesh,icolor)

  USE data_types
  USE GLOBAL_DATA 
  USE MYMPI

  IMPLICIT NONE

!---------------------------------------------------------------------------------------------------
! Dummy variables:

  TYPE(mesh_pointers),INTENT(IN) :: current_mesh
  INTEGER,INTENT(IN) :: icolor


! Local Variables

  REAL(KIND=double), DIMENSION(:,:,:), POINTER :: px

  INTEGER :: i,j,k,is,ks, snbr,rnbr, indx, tmpdnbr(4)
  INTEGER :: xe,ze,ye, elcount, ystride
  INTEGER :: mn, indsend(4),indrecv(4), nsend,nrecv
  INTEGER :: status(MPI_STATUS_SIZE),tmpcomm
!---------------------------------------------------------------------------------------------------

  px =>  current_mesh%x ; 
  xe = current_mesh%xe ; ye = current_mesh%ye; ze = current_mesh%ze;
  tmpdnbr(1:4) = current_mesh%dnbr(1:4)
  tmpcomm = current_mesh%comm3d
  mn = current_mesh%mesh_num
  indsend(1)=0;indsend(2)=xe;indsend(3)=0;indsend(4)=ze;
  indrecv(1)=-1;indrecv(2)=xe+1;indrecv(3)=-1;indrecv(4)=ze+1;
  ystride=ye+1

!
!   In rocfire domain decomposition, each processor has four neighbors
!   Consider here 2 at the time, first x after z
!
  IF(isparallel(1)) THEN      !parallel in x

     DO snbr =2,1,-1
        rnbr=rmap(snbr)
        nsend=current_mesh%nsnd(snbr,icolor)
        nrecv=current_mesh%nrcv(rnbr,icolor)

        i=indsend(snbr)
        do indx = 1,nsend
           k=current_mesh%sndclr(indx,snbr,icolor)
           do j = 0,ye 
              elcount=(indx-1)*ystride+j
              slambda(elcount) = px(i,k,j)
           enddo
        enddo


        CALL MPI_SENDRECV(slambda,nsend*ystride,&
             MPI_DOUBLE_PRECISION,tmpdnbr(snbr),0,rlambda,&
             nrecv*ystride,MPI_DOUBLE_PRECISION,tmpdnbr(rnbr),0,&
             tmpcomm,status,ierr)

        i=indrecv(rnbr)
        do indx = 1,nrecv
           k=current_mesh%rcvclr(indx,rnbr,icolor)
           do j = 0,ye   
              elcount=(indx-1)*ystride+j
              px(i,k,j) = rlambda(elcount)
           enddo
        enddo
     ENDDO

  ELSE

!
!
     DO snbr =1,2
        rnbr=rmap(snbr)
        nrecv=current_mesh%nrcv(rnbr,icolor)
        i=indrecv(rnbr);is=indsend(snbr)
        do indx = 1,nrecv
           k=current_mesh%rcvclr(indx,rnbr,icolor)
           do j = 0,ye
              px(i,k,j)=px(is,k,j)  
           enddo
        enddo
     ENDDO

  ENDIF

  IF(isparallel(2)) THEN      !parallel in z

     DO snbr =3,4
        rnbr=rmap(snbr)
        nsend=current_mesh%nsnd(snbr,icolor)
        nrecv=current_mesh%nrcv(rnbr,icolor)

        k=indsend(snbr)
        do indx = 1,nsend
           i=current_mesh%sndclr(indx,snbr,icolor)
           do j = 0,ye   
              elcount=(indx-1)*ystride+j
              slambda(elcount) = px(i,k,j)
           enddo
        enddo

        CALL MPI_SENDRECV(slambda,nsend*ystride,&
             MPI_DOUBLE_PRECISION,tmpdnbr(snbr),0,rlambda,&
             nrecv*ystride,MPI_DOUBLE_PRECISION,tmpdnbr(rnbr),0,&
             tmpcomm,status,ierr)

        k=indrecv(rnbr)
        do indx = 1,nrecv
           i=current_mesh%rcvclr(indx,rnbr,icolor)
           do j = 0,ye   
              elcount=(indx-1)*ystride+j
              px(i,k,j) = rlambda(elcount)
           enddo
        enddo
     ENDDO

  ELSE

     DO snbr =3,4
        rnbr=rmap(snbr)
        nrecv=current_mesh%nrcv(rnbr,icolor)
        k=indrecv(rnbr);ks=indsend(snbr)
        do indx = 1,nrecv
           i=current_mesh%rcvclr(indx,rnbr,icolor)
           do j = 0,ye
              px(i,k,j)=px(i,ks,j)  
           enddo
        enddo
     ENDDO

  ENDIF

!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!.................................................................................
  IF(issymmetric) THEN
     IF(dSerialRange(1).eq.0)         px(-1,:,:) =   px(1,:,:)
     IF(coords(1)+1 == pdims(1)) px(xe+1,:,:) = px(xe,:,:)
     IF(is2D) THEN
        px(:,-1,:) = px(:,0,:)
        px(:,-2,:) = px(:,0,:)
        px(:,+1,:) = px(:,0,:)
        px(:,+2,:) = px(:,0,:)
     ELSE
        IF(dSerialRange(3).eq.0)         px(:,-1,:) =   px(:,1,:)
        IF(coords(2)+1 == pdims(2)) px(:,ze+1,:) = px(:,ze,:)
     ENDIF

  ENDIF

!----------------------------------------------------------------------------------
  RETURN
END SUBROUTINE color_swap
!**********************************************************************************


!**********************************************************************************
SUBROUTINE fill_voids(current_mesh,fine_mesh)

  USE data_types
  USE GLOBAL_DATA 
  USE MYMPI

  IMPLICIT NONE

!---------------------------------------------------------------------------------------------------
! Dummy variables:

  TYPE(mesh_pointers),INTENT(IN) :: current_mesh,fine_mesh

!   LOCAL variables

  REAL(KIND=double), DIMENSION(:,:,:), POINTER :: px

  INTEGER :: i,j,k,tmpdnbr(4),tmpcomm
  INTEGER :: nxc,nzc,nyc,nxf,nzf
  INTEGER :: xs,zs,xe,ze,ye, elcount
  INTEGER :: status(MPI_STATUS_SIZE)
!---------------------------------------------------------------------------------------------------


  px =>  current_mesh%x ; 
  nxc = current_mesh%xe+1 ; nyc = current_mesh%ye+1; nzc = current_mesh%ze+1;
  nxf = fine_mesh%xe+1 ;nzf = fine_mesh%ze+1;
  tmpdnbr(1:4) = fine_mesh%dnbr(1:4)
  tmpcomm = fine_mesh%comm3d

  do i = 1,2
     if(current_mesh%recv_void(i) == 1) then
        elcount = 2*(nzc+2)*nyc
        call MPI_RECV(rlambda,elcount,MPI_DOUBLE_PRECISION,tmpdnbr(i),0,tmpcomm,status,ierr)
        elcount = 0
        do j=0,nyc-1
           do k=-1,nzc
              px(-1,k,j) =  rlambda(elcount)
              px(nxc,k,j) = rlambda(elcount+1)
              elcount = elcount+2
           enddo
        enddo
     endif
  enddo
!
  do k = 3,4
     if(current_mesh%recv_void(k) == 1) then
        elcount = 2*(nxc+2)*nyc
        call MPI_RECV(rlambda,elcount,MPI_DOUBLE_PRECISION,tmpdnbr(k),0,tmpcomm,status,ierr)
        elcount = 0
        do j=0,nyc-1
           do i=-1,nxc
              px(i,-1,j) =  rlambda(elcount)
              px(i,nzc,j) = rlambda(elcount+1)
              elcount = elcount+2
           enddo
        enddo
     endif
  enddo

!.........................................................................................

  if(current_mesh%send_void(3) == 1) then
     elcount = 0;zs=-1
     do j=0,nyc-1
        do i=-1,nxc
           slambda(elcount) = px(i,zs,j)
           slambda(elcount+1) = px(i,zs+1,j)
           elcount = elcount+2
        enddo
     enddo
     call MPI_SEND(slambda,elcount,MPI_DOUBLE_PRECISION,tmpdnbr(3),0,tmpcomm,status,ierr)
  endif

  if(current_mesh%send_void(4) == 1 ) then
     elcount = 0;zs=nzc-1
     do j=0,nyc-1
        do i=-1,nxc
           slambda(elcount) = px(i,zs,j)
           slambda(elcount+1) = px(i,zs+1,j)
           elcount = elcount+2
        enddo
     enddo
     call MPI_SEND(slambda,elcount,MPI_DOUBLE_PRECISION,tmpdnbr(4),0,tmpcomm,status,ierr)
  endif
!
  if(current_mesh%send_void(1) == 1 ) then
     elcount = 0;xs=-1
     do j=0,nyc-1
        do k=-1,nzc
           slambda(elcount) = px(xs,k,j)
           slambda(elcount+1) = px(xs+1,k,j)
           elcount = elcount+2
        enddo
     enddo
     call MPI_SEND(slambda,elcount,MPI_DOUBLE_PRECISION,tmpdnbr(1),0,tmpcomm,status,ierr)
  endif

  if(current_mesh%send_void(2) == 1 ) then
     elcount = 0;xs=nxc-1
     do j=0,nyc-1
        do k=-1,nzc
           slambda(elcount) = px(xs,k,j)
           slambda(elcount+1) = px(xs+1,k,j)
           elcount = elcount+2
        enddo
     enddo
     call MPI_SEND(slambda,elcount,MPI_DOUBLE_PRECISION,tmpdnbr(2),0,tmpcomm,status,ierr)
  endif
!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!.................................................................................
  xe = current_mesh%xe
  ze = current_mesh%ze
  IF(issymmetric) THEN
     IF(dSerialRange(1).eq.0)         px(-1,:,:) =   px(1,:,:)
     IF(coords(1)+1 == pdims(1)) px(xe+1,:,:) = px(xe,:,:)
     IF(dSerialRange(3).eq.0)         px(:,-1,:) =   px(:,1,:)
     IF(coords(2)+1 == pdims(2)) px(:,ze+1,:) = px(:,ze,:)
  ENDIF


!----------------------------------------------------------------------------------
  RETURN
END SUBROUTINE fill_voids
!**********************************************************************************

!**********************************************************************************
SUBROUTINE fill_voids4(current_mesh,fine_mesh)

  USE data_types
  USE GLOBAL_DATA 
  USE MYMPI

  IMPLICIT NONE

!---------------------------------------------------------------------------------------------------
! Dummy variables:

  TYPE(mesh_pointers),INTENT(IN) :: current_mesh,fine_mesh

!   LOCAL variables

  REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: px

  INTEGER :: i,j,k,tmpdnbr(4),tmpcomm
  INTEGER :: eqn1,neqn,eqn
  INTEGER :: nxc,nzc,nyc,nxf,nzf,ystride
  INTEGER :: xs,zs,xe,ze,ye, elcount
  INTEGER :: status(MPI_STATUS_SIZE)
!---------------------------------------------------------------------------------------------------


  px =>  current_mesh%w ; 
  nxc = current_mesh%xe+1 ; nyc = current_mesh%ye+1; nzc = current_mesh%ze+1;
  nxf = fine_mesh%xe+1 ;nzf = fine_mesh%ze+1;
  tmpdnbr(1:4) = fine_mesh%dnbr(1:4)
  tmpcomm = fine_mesh%comm3d
  eqn1=1; neqn=neqmax-1
  ystride =nyc+1*(neqn - eqn1+1)

  do i = 1,2
     if(current_mesh%recv_void(i) == 1) then
        elcount = 2*(nzc+2)*ystride
        call MPI_RECV(rlambda,elcount,MPI_DOUBLE_PRECISION,tmpdnbr(i),0,tmpcomm,status,ierr)
        elcount = 0
        do j=0,nyc
           do k=-1,nzc
              do eqn = eqn1,neqn
                 px(-1,k,j,eqn) =  rlambda(elcount)
                 px(nxc,k,j,eqn) = rlambda(elcount+1)
                 elcount = elcount+2
              enddo
           enddo
        enddo
     endif
  enddo
!
  do k = 3,4
     if(current_mesh%recv_void(k) == 1) then
        elcount = 2*(nxc+2)*ystride
        call MPI_RECV(rlambda,elcount,MPI_DOUBLE_PRECISION,tmpdnbr(k),0,tmpcomm,status,ierr)
        elcount = 0
        do j=0,nyc
           do i=-1,nxc
              do eqn = eqn1,neqn
                 px(i,-1,j,eqn) =  rlambda(elcount)
                 px(i,nzc,j,eqn) = rlambda(elcount+1)
                 elcount = elcount+2
              enddo
           enddo
        enddo
     endif
  enddo

!.........................................................................................

  if(current_mesh%send_void(3) == 1) then
     elcount = 0;zs=-1
     do j=0,nyc
        do i=-1,nxc
           do eqn = eqn1,neqn
              slambda(elcount) = px(i,zs,j,eqn)
              slambda(elcount+1) = px(i,zs+1,j,eqn)
              elcount = elcount+2
           enddo
        enddo
     enddo
     call MPI_SEND(slambda,elcount,MPI_DOUBLE_PRECISION,tmpdnbr(3),0,tmpcomm,status,ierr)
  endif

  if(current_mesh%send_void(4) == 1 ) then
     elcount = 0;zs=nzc-1
     do j=0,nyc
        do i=-1,nxc
           do eqn = eqn1,neqn
              slambda(elcount) = px(i,zs,j,eqn)
              slambda(elcount+1) = px(i,zs+1,j,eqn)
              elcount = elcount+2
           enddo
        enddo
     enddo
     call MPI_SEND(slambda,elcount,MPI_DOUBLE_PRECISION,tmpdnbr(4),0,tmpcomm,status,ierr)
  endif
!
  if(current_mesh%send_void(1) == 1 ) then
     elcount = 0;xs=-1
     do j=0,nyc
        do k=-1,nzc
           do eqn = eqn1,neqn
              slambda(elcount) = px(xs,k,j,eqn)
              slambda(elcount+1) = px(xs+1,k,j,eqn)
              elcount = elcount+2
           enddo
        enddo
     enddo
     call MPI_SEND(slambda,elcount,MPI_DOUBLE_PRECISION,tmpdnbr(1),0,tmpcomm,status,ierr)
  endif

  if(current_mesh%send_void(2) == 1 ) then
     elcount = 0;xs=nxc-1
     do j=0,nyc
        do k=-1,nzc
           do eqn = eqn1,neqn
              slambda(elcount) = px(xs,k,j,eqn)
              slambda(elcount+1) = px(xs+1,k,j,eqn)
              elcount = elcount+2
           enddo
        enddo
     enddo
     call MPI_SEND(slambda,elcount,MPI_DOUBLE_PRECISION,tmpdnbr(2),0,tmpcomm,status,ierr)
  endif

!----------------------------------------------------------------------------------
  RETURN
END SUBROUTINE fill_voids4
!**********************************************************************************

!**********************************************************************************
SUBROUTINE timer(ttime,flag)

  USE data_types
  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!---------------------------------------------------------------------------------------------------
! Dummy variables:

  REAL *8,INTENT(INOUT) :: ttime(2)
  INTEGER ,INTENT(IN) :: flag         !0 start the timer 1 stops the timer

!   LOCAL variables
!---------------------------------------------------------------------------------------------------

  if (flag == 0 ) then
     ttime(2) = MPI_WTIME()
  else
     ttime(2) = MPI_WTIME() - ttime(2)
     ttime(1) = ttime(1) + ttime(2)
  endif
!----------------------------------------------------------------------------------
  RETURN
END SUBROUTINE timer
!
!*************************************************************************
SUBROUTINE twolayer_swap(col,current_mesh)

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE    

!------------------------------------------------------------------------
! Dummy variables:
  TYPE(mesh_pointers),INTENT(IN) :: current_mesh
  INTEGER, INTENT(IN) :: col(4)

! Local variables
  REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: pw

  INTEGER ::  i,k,j,xs,xe,ys,ye,zs,ze, r1,r2,layer, s1,s2,iss(maxsize)
  INTEGER ::  status(MPI_STATUS_SIZE),eqn,eqn1,neqn,maxmem
  INTEGER ::  indsend(4,2),indrecv(4,2), snbr,rnbr, elcount, elrecv
  INTEGER ::  tmpdnbr(4),tmpcomm

  REAL*8  ::  ss(maxsize)
!--------------------------------------------------------------------------

!
! THIS IS the array being swapped
!
  pw =>  current_mesh%w; 
!
  eqn1 = col (1)
  neqn = col (2)
  ys = col(3)
  ye = col(4)
  xe = current_mesh%xe ; ze = current_mesh%ze;
  tmpdnbr(1:4) = current_mesh%dnbr(1:4)
  tmpcomm = current_mesh%comm3d

  indsend(1,1)=0;indsend(1,2) = 1; indsend(2,1)=xe-1;indsend(2,2) = xe;
  indsend(3,1)=0;indsend(3,2) = 1; indsend(4,1) = ze -1;indsend(4,2)=ze;

  indrecv(1,1) = -2;indrecv(1,2)=-1; indrecv(2,1)=xe+1; indrecv(2,2) = xe+2;
  indrecv(3,1) = -2;indrecv(3,2)=-1; indrecv(4,1)=ze+1; indrecv(4,2) = ze+2;

10 continue

  maxmem=min(minval(shape(slambda)),minval(shape(rlambda)))


  IF(isparallel(1)) THEN      !parallel in x

     DO snbr =2,1,-1
        rnbr=rmap(snbr)

        s1=indsend(snbr,1); s2 = indsend(snbr,2); layer = s2-s1
        elcount = 0
        DO k = 0, ze
           do j = ys,ye
              do eqn = eqn1,neqn
                 slambda(elcount:elcount+layer) = pw(s1:s2,k,j,eqn)
                 elcount=elcount+layer+1
              enddo
           enddo
        enddo
!
! check the dimension of the buffer used in the communication
!
        if(elcount.gt.maxmem) THEN 
           write(*,*) 'INSUFFICIENT MEMORY IN COMMUNICATION',elcount,maxmem
           DEALLOCATE(slambda,rlambda)
           ALLOCATE(slambda(0:elcount+1),rlambda(0:elcount+1))
           GOTO 10
        endif

        elrecv=elcount
        CALL MPI_SENDRECV(slambda,elcount,&
             MPI_DOUBLE_PRECISION,dnbr(snbr),0,rlambda,&
             elrecv,MPI_DOUBLE_PRECISION,dnbr(rnbr),0,&
             tmpcomm,status,ierr)

        r1=indrecv(rnbr,1);r2 =indrecv(rnbr,2); layer = r2-r1
        elcount=0
        DO k = 0, ze
           do j = ys,ye
              do eqn = eqn1,neqn
                 pw(r1:r2,k,j,eqn) = rlambda(elcount:elcount+layer)
                 elcount=elcount+layer+1
              enddo
           enddo
        enddo

     ENDDO

  ELSE
!
     DO snbr = 1,2
        rnbr=rmap(snbr)
        r1=indrecv(rnbr,1);r2 = indrecv(rnbr,2);
        s1=indsend(snbr,1);s2 = indsend(snbr,2);
        DO k = 0, ze
           do j = ys,ye
              do eqn = eqn1,neqn
                 pw(r1:r2,k,j,eqn)=pw(s1:s2,k,j,eqn)
              enddo
           enddo
        enddo
     ENDDO

  ENDIF

  IF(isparallel(2)) THEN   !parallel in z

     DO snbr =4,3,-1
        rnbr=rmap(snbr)

        s1=indsend(snbr,1); s2 = indsend(snbr,2); layer = s2-s1
        elcount = 0
        DO i = -2, xe+2
           do j = ys,ye
              do eqn = eqn1,neqn
                 slambda(elcount:elcount+layer) = pw(i,s1:s2,j,eqn)
                 elcount=elcount+layer+1
              enddo
           enddo
        enddo
!
! check the dimension of the buffer used in the communication
!
        if(elcount.gt.maxmem) THEN 
           write(*,*) 'INSUFFICIENT MEMORY IN COMMUNICATION',elcount,maxmem
           DEALLOCATE(slambda,rlambda)
           ALLOCATE(slambda(0:elcount+1),rlambda(0:elcount+1))
           GOTO 10
        endif

        elrecv=elcount
        CALL MPI_SENDRECV(slambda,elcount,&
             MPI_DOUBLE_PRECISION,dnbr(snbr),0,rlambda,&
             elrecv,MPI_DOUBLE_PRECISION,dnbr(rnbr),0,&
             tmpcomm,status,ierr)

        r1=indrecv(rnbr,1);r2 =indrecv(rnbr,2); layer = r2-r1
        elcount=0
        do i = -2, xe+2
           do j = ys,ye
              do eqn = eqn1,neqn
                 pw(i,r1:r2,j,eqn) = rlambda(elcount:elcount+layer)
                 elcount=elcount+layer+1
              enddo
           enddo
        enddo

     ENDDO

  ELSE
!
     DO snbr = 3,4
        rnbr=rmap(snbr)
        r1=indrecv(rnbr,1);r2 = indrecv(rnbr,2);
        s1=indsend(snbr,1);s2 = indsend(snbr,2);
        do i = -2, xe+2
           do j = ys,ye
              do eqn = eqn1,neqn
                 pw(i,r1:r2,j,eqn)=pw(i,s1:s2,j,eqn)
              enddo
           enddo
        enddo
     ENDDO
!
  ENDIF

!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!.................................................................................
  !!indsend(1,1)=0;indsend(1,2) = 1; indsend(2,1)=xe-1;indsend(2,2) = xe;
  !!indsend(3,1)=0;indsend(3,2) = 1; indsend(4,1) = ze -1;indsend(4,2)=ze;
  !!
  !!indrecv(1,1) = -2;indrecv(1,2)=-1; indrecv(2,1)=xe+1; indrecv(2,2) = xe+2;
  !!indrecv(3,1) = -2;indrecv(3,2)=-1; indrecv(4,1)=ze+1; indrecv(4,2) = ze+2;

  xe = current_mesh%xe
  ze = current_mesh%ze


  IF(issymmetric) THEN
     DO eqn = eqn1,neqn


        ss = one
        iss = 0
        IF(neqn <= 3) THEN 
           ss(1) = -1.0d0
           iss(1) = -1
        ENDIF

        IF(dSerialRange(1).eq.0)           pw(-2,:,:,eqn) =   ss(eqn)*pw(2+iss(eqn),:,:,eqn)
        IF(dSerialRange(1).eq.0)           pw(-1,:,:,eqn) =   ss(eqn)*pw(1+iss(eqn),:,:,eqn)
        IF(coords(1)+1 == pdims(1) )  pw(xe+1,:,:,eqn) = ss(eqn)*pw(xe+iss(eqn),:,:,eqn)
        IF(coords(1)+1 == pdims(1) )  pw(xe+2,:,:,eqn) = ss(eqn)*pw(xe-1+iss(eqn),:,:,eqn)

        IF(is2D) THEN

           IF(dSerialRange(3).eq.0)           pw(:,-1,:,eqn) =   pw(:,1,:,eqn)
           IF(dSerialRange(3).eq.0)           pw(:,-2,:,eqn) =   pw(:,2,:,eqn)
           IF(coords(2)+1 == pdims(2) )  pw(:,ze+1,:,eqn) = pw(:,ze,:,eqn)
           IF(coords(2)+1 == pdims(2) )  pw(:,ze+2,:,eqn) = pw(:,ze-1,:,eqn)

        ELSE


           IF(dSerialRange(3).eq.0)           pw(:,-2,:,eqn) =   ss(eqn)*pw(:,2+iss(eqn),:,eqn)
           IF(dSerialRange(3).eq.0)           pw(:,-1,:,eqn) =   ss(eqn)*pw(:,1+iss(eqn),:,eqn)
           IF(coords(2)+1 == pdims(2) )  pw(:,ze+1,:,eqn) = ss(eqn)*pw(:,ze+iss(eqn),:,eqn)
           IF(coords(2)+1 == pdims(2) )  pw(:,ze+2,:,eqn) = ss(eqn)*pw(:,ze-1+iss(eqn),:,eqn)
        ENDIF

!
!   Only for 2D: MAKE a Z swap just to make sure
!
        if(ipack == 0) then
           IF(dserialrange(3).eq.0)           pw(:,-1,:,eqn) =   pw(:,0,:,eqn)
           IF(dserialrange(3).eq.0)           pw(:,-2,:,eqn) =   pw(:,0,:,eqn)
           IF(coords(2)+1 == pdims(2) )  pw(:,ze+1,:,eqn) = pw(:,ze,:,eqn)
           IF(coords(2)+1 == pdims(2) )  pw(:,ze+2,:,eqn) = pw(:,ze,:,eqn)
        endif

     ENDDO
  ENDIF


!------------------------------------------------------------------------
  RETURN
END SUBROUTINE twolayer_swap
!**********************************************************************************
