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
  
!  IF(ddrange(1).eq.0) xs = 1
!  IF(ddrange(2).eq.nx) xe = xe -1
!  IF(ddrange(3).eq.0) zs = 1
!  IF(ddrange(4).eq.nz) ze = ze - 1
  

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
  xe = drange(2)
  ze = drange(4)
  IF(issymmetric) THEN
   IF(ddrange(1).eq.0)              psi(-1,:,:) =     psi(1,:,:)
   IF(coords(1)+1 == pdims(1) )     psi(xe+1,:,:) =   psi(xe,:,:)
   IF(ddrange(3).eq.0)              psi(:,-1,:) =     psi(:,1,:)
   IF(coords(2)+1 == pdims(2) )     psi(:,ze+1,:) =   psi(:,ze,:)
  ENDIF
  
!-----------------------------------------------------
  RETURN
END SUBROUTINE FILL_PSI_GHOST
!********************************************************************

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
!     IF(ddrange(1).eq.0) xsend(1) = 1
!     IF(ddrange(2).eq.nx) xsend(2) = xe -1
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
!     IF(ddrange(3).eq.0) zsend(1) = 1
!     IF(ddrange(4).eq.nz) zsend(2) = ze - 1
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
  xe = drange(2)
  ze = drange(4)
  IF(issymmetric) THEN
   IF(ddrange(1).eq.0)      phi(-1,:) =     phi(1,:)
   IF(ddrange(1).eq.0)      phi(-2,:) =     phi(2,:)
   IF(coords(1)+1 == pdims(1) ) THEN
      phi(xe+1,:) =   phi(xe,:)
      phi(xe+2,:) =   phi(xe-1,:)
   ENDIF
   IF(ddrange(3).eq.0)      phi(:,-1) =     phi(:,1)
   IF(ddrange(3).eq.0)      phi(:,-2) =     phi(:,2)
   IF(coords(2)+1 == pdims(2) ) THEN
      phi(:,ze+1) =   phi(:,ze)
      phi(:,ze+2) =   phi(:,ze-1)
   ENDIF
  ENDIF
  
  RETURN
END SUBROUTINE FILL_PHI_GHOST


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
  LOGICAL :: do_not_receive(2)
!------------------------------------------------------------------------------------------------

  xs = 0
  xe = col(2)
  zs = 0
  ze = col(4)
  iadd=0;kadd=0;
  tmpdnbr(1:4) = current_mesh%dnbr(1:4)
  tmpcomm = current_mesh%comm3d
 

  do_not_receive(1) = is_proc_bnd(1) .AND. .NOT. is_periodic_xz(0)
  do_not_receive(2) = is_proc_bnd(2) .AND. .NOT. is_periodic_xz(0)
  do_not_receive = .FALSE.

  IF(isparallel(1)) THEN      !parallel in x

     elcount = 0
     DO k = col(3), col(4)
        DO j = col(5), col(6)
           slambda(elcount) = px(xe,k,j)    !SEND the last colon, kadd
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
           if(.NOT. do_not_receive(1)) &
                px(xs-1,k,j) = rlambda(elcount)   !RECEIVE the last colon, kadd
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
           if(.NOT. do_not_receive(2)) &
                px(xe+1,k,j) = rlambda(elcount)    !RECEIVE the first colon
           elcount = elcount + 1
        END DO
     END DO


  ELSEIF(.NOT. ANY(do_not_receive)) THEN
     DO k = col(3), col(4)
     DO j = col(5), col(6)
        px(xs-1,k,j) = px(xe,k,j)
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
     write(*,*) mywid,':: Doing a Y send?!?'
  END IF


  ye = current_mesh%ye
  ys = 0
  if(is_periodic_y) then
     px(:,:,ye+1) = px(:,:,ys)
     px(:,:,ys-1) = px(:,:,ye)
  endif

!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!.................................................................................  

  IF(issymmetric) THEN
   IF(ddrange(1).eq.0)         px(-1,:,:) =   px(1,:,:)
   IF(coords(1)+1 == pdims(1)) px(xe+1,:,:) = px(xe,:,:)
   IF(ddrange(3).eq.0)         px(:,-1,:) =   px(:,1,:)
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

  INTEGER ::  i,k,j,xs,xe,ys,ye,zs,ze,icolor,color_start,color_end, indx
  INTEGER ::  status(MPI_STATUS_SIZE),xrange,yrange,size,eqn,eqn1,neqn,maxmem
  INTEGER ::  indsend(4),indrecv(4), nsend,nrecv, snbr,rnbr, is,ks, elcount, elrecv
  INTEGER ::  tmpdnbr(4),tmpcomm, iss(neqmax)
  REAL(KIND=double) ::  ss(neqmax)
  LOGICAL :: do_not_receive(2)
!--------------------------------------------------------------------------

  if (col(1) > 0) then
     pw =>  current_mesh%w; 
  else
     pw => wrk4
  endif
!
  eqn_ref = col(7)
!
  color_start= abs(col(1))
  color_end= col(2)
  eqn1 = col (3)
  neqn = col (4)
  ys = col(5)
  ye = col(6)
  xe = current_mesh%xe; 
  ze = current_mesh%ze;



  tmpdnbr(1:4) = current_mesh%dnbr(1:4)
  tmpcomm = current_mesh%comm3d
  indsend(1)=0;indsend(2)=xe;indsend(3)=0;indsend(4)=ze;
  indrecv(1)=-1;indrecv(2)=xe+1;indrecv(3)=-1;indrecv(4)=ze+1;

  do_not_receive(1) = is_proc_bnd(1) .AND. .NOT. is_periodic_xz(eqn_ref)
  do_not_receive(2) = is_proc_bnd(2) .AND. .NOT. is_periodic_xz(eqn_ref)

10 continue

  maxmem=min(minval(shape(slambda)),minval(shape(rlambda)))


  IF(isparallel(1)) THEN      !parallel in x

    DO snbr =2,1,-1
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

       if(do_not_receive(rnbr)) GOTO 100

       i=indrecv(rnbr)
       elcount=0
       do icolor = color_start,color_end 
         nrecv= current_mesh%nrcv(rnbr,icolor)
         do indx = 1,nrecv
            k= current_mesh%rcvclr(indx,rnbr,icolor)
            do j = ys,ye
            do eqn = eqn1,neqn
               pw(i,k,j,eqn) = rlambda(elcount)
               elcount=elcount+1
            enddo
            enddo
         enddo
       enddo

100 CONTINUE
  
    ENDDO


  ELSEIF(.NOT. ANY(do_not_receive) )THEN
!
    DO snbr =1,2
      rnbr=rmap(snbr)
      i=indrecv(rnbr);is=indsend(snbr)
      do icolor = color_start,color_end
         nrecv= current_mesh%nrcv(rnbr,icolor)
         do indx = 1,nrecv
            k= current_mesh%rcvclr(indx,rnbr,icolor)
            do j = ys,ye
            do eqn = eqn1,neqn
               pw(i,k,j,eqn)=pw(is,k,j,eqn)
            enddo
            enddo
         enddo
      enddo
    ENDDO

  ENDIF  
  
  IF(isparallel(2) .AND. .NOT. is2d) THEN   !parallel in z

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

       k=indrecv(rnbr)
       elcount=0
       do icolor = color_start,color_end 
         nrecv= current_mesh%nrcv(rnbr,icolor)
         do indx = 1,nrecv
            i= current_mesh%rcvclr(indx,rnbr,icolor)
            do j = ys,ye
            do eqn = eqn1,neqn
               pw(i,k,j,eqn) = rlambda(elcount)
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
      k=indrecv(rnbr);ks=indsend(snbr)
      do icolor = color_start,color_end
         nrecv= current_mesh%nrcv(rnbr,icolor)
         do indx = 1,nrecv
            i= current_mesh%rcvclr(indx,rnbr,icolor)
            do j = ys,ye
            do eqn = eqn1,neqn
               pw(i,k,j,eqn)=pw(i,ks,j,eqn)
            enddo
            enddo
         enddo
      enddo
    ENDDO

 ENDIF  

 if(is_periodic_y) then
    ye = current_mesh%ye
    ys = 0
    if(current_mesh%mesh_num == 1) then  
       ye = ny-1
       ys = 1
    endif
    pw(:,:,ye+1,:) = pw(:,:,ys,:)
    pw(:,:,ys-1,:) = pw(:,:,ye,:)
    pw(:,:,ye+2,:) = pw(:,:,ys+1,:)
    pw(:,:,ys-2,:) = pw(:,:,ye-1,:)
 endif
!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!.................................................................................  

  IF(issymmetric) THEN
   DO eqn = eqn1,neqn

!     ss=min(dble(2*eqn) - 3.0d0 , one )    !eqn=1,2,3  ss=-1,1,1   iss=-1,0,0
!     iss = min(2*eqn - 3 , 0)
     ss = (/ -1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0 /)
     iss = (/ -1 , 0 , 0 , 0 , 0, 0/)

     IF(var_updated == 0) THEN 
        ss = one
        iss = 0
     ENDIF

     IF(ddrange(1).eq.0)           pw(-1,:,:,eqn) =   ss(eqn) * pw(1+iss(eqn),:,:,eqn)
     IF(coords(1)+1 == pdims(1) )  pw(xe+1,:,:,eqn) = ss(eqn) * pw(xe+iss(eqn),:,:,eqn)

!---2D only

     pw(:,-1,:,eqn) =   pw(:,0,:,eqn)
     pw(:,ze+1,:,eqn) = pw(:,ze,:,eqn)
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

    TYPE(mesh_pointers),INTENT(INOUT) :: current_mesh
    INTEGER,INTENT(IN) :: icolor
    
   
  ! Local Variables

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: px

    INTEGER :: iii,ii,i,j,k,is,ks, snbr,rnbr, indx, tmpdnbr(4)
    INTEGER :: xe,ze,ye,ys, elcount, ystride,itemp
    INTEGER :: mn, indsend(4),indrecv(4), nsend,nrecv
    INTEGER :: status(MPI_STATUS_SIZE),tmpcomm
    LOGICAL :: do_not_receive(2)
!---------------------------------------------------------------------------------------------------

    px =>  current_mesh%x ; 
    xe = current_mesh%xe ; ye = current_mesh%ye; ze = current_mesh%ze;
    tmpdnbr(1:4) = current_mesh%dnbr(1:4)
    tmpcomm = current_mesh%comm3d
    mn = current_mesh%mesh_num
    indsend(1)=0;indsend(2)=xe;indsend(3)=0;indsend(4)=ze;
    indrecv(1)=-1;indrecv(2)=xe+1;indrecv(3)=-1;indrecv(4)=ze+1;
    ystride=ye+1

    do_not_receive(1) = is_proc_bnd(1) .AND. .NOT. is_periodic_xz(4)
    do_not_receive(2) = is_proc_bnd(2) .AND. .NOT. is_periodic_xz(4)

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

         if(do_not_receive(rnbr)) CYCLE
         
         i=indrecv(rnbr)
         do indx = 1,nrecv
            k=current_mesh%rcvclr(indx,rnbr,icolor)
         do j = 0,ye   
            elcount=(indx-1)*ystride+j
            px(i,k,j) = rlambda(elcount)
         enddo
         enddo
      ENDDO
   
    ELSEIF(.NOT. ANY(do_not_receive) )THEN

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

    ye = current_mesh%ye
    ys = 0
    if(is_periodic_y) then
       px(:,:,ye+1) = px(:,:,ys)
       px(:,:,ys-1) = px(:,:,ye)
    endif

  call DIRICHLET_BC(itemp,current_mesh)
!.................................................................................
! MODIFICATIONS FOR SYMMETRIC BOUNDARIES
!................................................................................. 

  IF(issymmetric) THEN
    IF(indsend(1)+ib_mg(mn) == 0)            px(-1,:,:)   = px(+1,:,:)
    IF(indsend(2)+ib_mg(mn) == nx_tot(mn)-1) px(xe+1,:,:) = px(xe,:,:)
    IF(is2D) THEN
      px(:,-1,:) = px(:,0,:)
      px(:,-2,:) = px(:,0,:)
      px(:,+1,:) = px(:,0,:)
      px(:,+2,:) = px(:,0,:)
    ELSE
       IF(indsend(3)+kb_mg(mn) == 0)            px(:,-1,:) =   px(:,1,:)
       IF(indsend(4)+kb_mg(mn) == nz_tot(mn)-1) px(:,ze+1,:) = px(:,ze,:)
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
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
   IF(issymmetric) THEN
     IF(ddrange(1).eq.0)         px(-1,:,:) =   px(1,:,:)
     IF(coords(1)+1 == pdims(1)) px(xe+1,:,:) = px(xe,:,:)
     IF(ddrange(3).eq.0)         px(:,-1,:) =   px(:,1,:)
     IF(coords(2)+1 == pdims(2)) px(:,ze+1,:) = px(:,ze,:)
   ENDIF

      
!----------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE fill_voids
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
!*************************************************************************
!

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
!     IF(ddrange(1).eq.0) xs = 1
!     IF(ddrange(2).eq.nx) xe = xe - 1
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
!     IF(ddrange(3).eq.0) zs = 1
!     IF(ddrange(4).eq.nz) ze = ze - 1
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
  xe = drange(2)
  ze = drange(4)
  IF(issymmetric) THEN
   IF(ddrange(1).eq.0)          f(-1,:,:,eqn) =   f(1,:,:,eqn)
   IF(coords(1)+1 == pdims(1) ) f(xe+1,:,:,eqn) = f(xe,:,:,eqn)
   IF(ddrange(3).eq.0)          f(:,-1,:,eqn) =   f(:,1,:,eqn)
   IF(coords(2)+1 == pdims(2) ) f(:,ze+1,:,eqn) = f(:,ze,:,eqn)
  ENDIF
  

  RETURN
END SUBROUTINE FILL_F_GHOST

!
!
! ********************************************************************
!
