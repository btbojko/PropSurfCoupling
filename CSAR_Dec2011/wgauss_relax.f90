!     ****************************************************************
!     *                                                              *
!     *                      subroutine gauss_relax                  *
!     *                                                              *
!     *                             solve Ax=f                       *
!     *                                                              *
!     *     with red and black decomposition in the x &  z           *
!     *                                                              *
!     *                                                              * 
!     ****************************************************************
!***************************************************************************************************
!NOTE:

!***************************************************************************************************

SUBROUTINE gauss_relax4 (current_mesh, nu, resid_flag)

  USE data_types
  USE GLOBAL_DATA 
  USE LUSOLVER, ONLY : LUSOLVE, LUINV
  USE MYMPI

  IMPLICIT NONE

  include "parallel.h"

!---------------------------------------------------------------------------------------------------
!   Dummy variables:

    TYPE(mesh_pointers) :: current_mesh

    INTEGER, INTENT(IN) :: nu                ! number of relaxations 
    LOGICAL, INTENT(IN) :: resid_flag        ! if this is true it evaluates the residual
   
!   Local Variables
    REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: px
    REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: rhs


    INTEGER :: icount,i,j,k,xe,ze,ys,ye,m,n,col(7),jj
    INTEGER :: indx,icolor,rcolor, tmpcomm
    INTEGER :: eqn1,neqn,eqn,mn,jp1,iteqn,js,jmax,n_samples
    REAL (KIND=double) :: coe,c2,c3,c4,c5,c8,Ts
    REAL (KIND=double) :: f1,f2,f3,f4,f5,gxy
    REAL (KIND=double) :: deltaf,addrate,rr,rhsg,myres(2),resmax
    REAL (KIND=double) :: c11,sumc,omg_txy,one_m_omg, sum(2), mysum(2)
    REAL (KIND=double) :: ftmp(neqgas), OLD_TEMP
!---------------------------------------------------------------------------------------------------

!
!   RETURN immediatly if the processor has no points
!
    if(current_mesh%blank) RETURN

    px => current_mesh%w ; rhs => current_mesh%r

    eqn_ref = 0
    xe = nxv(eqn_ref) ; ze = current_mesh%ze
    ys = 1; ye = nyv(eqn_ref)

    eqn1 = 1; neqn=neqmax
    col(5)=0;col(6)=ye;
    col(7) = eqn_ref
    mn = current_mesh%mesh_num

    tmpcomm = current_mesh%comm3d
    omg_txy=1.00d0
    one_m_omg = (one - omg_txy)

!--- Start and end of line relaxation loops
    jmax = ny
    js = 0 
!
!   START ITERATIONS
! 
    col(3)=eqn1;col(4)=neqn;
    do 25  icount = 1,nu


    do 24 icolor = 1,2


       rcolor=mod(icolor+2,2)+1     !2,1
       col(1:2) = rcolor
       call parallel_swap4(col,current_mesh)


       DO indx = 1 , eqn_nclr(icolor,0)

          i = eqn_iclr(indx,icolor,0)
          k = eqn_kclr(indx,icolor,0)

          if(.NOT. saved_LU_vars) then
             do m = 1,neqgas
                do n = 1,neqgas
                   do j = js,jmax
                      Bb(m,n,j) = Bm(m,n,i,k,j)
                   enddo
                enddo
             enddo
          endif

          fyv(1:neqgas,js) = px(i,k,js,1:neqgas)
          ay(js) = 0.0d0
          cy(js) = 0.0d0
          do j=ys,ye
             if(j<ye)then
                jp1 = j+1
             else
                jp1 = j
             endif
             c2=coeff(2,i,k,j)
             c3=coeff(3,i,k,j)
             c8=coeff(8,i,k,j)
             ay(j) = - coeff(4,i,k,j)
             cy(j) = - coeff(5,i,k,j)

             do eqn = 1,neqgas
                f2=px(i-1,k,j,eqn)
                f3=px(i+1,k,j,eqn)
                gxy = (px(i+1,k,jp1,eqn)-px(i+1,k,j-1,eqn)-px(i-1,k,jp1,eqn)+px(i-1,k,j-1,eqn))
                fyv(eqn,j) = dfdt(i,k,j,eqn) + c2*f2+c3*f3 + c8*gxy
             enddo   !eqns

          enddo  !j
          fyv(1:neqgas,jmax) = 0.0d0
          ay(jmax) = -1.0d0
          cy(jmax) = 0.0d0

! --- THOMAS algorithm--start
          do j=js+1,jmax-1
             if(.NOT. saved_LU_vars) then
                CALL LUINV(Bb(1:neqgas,1:neqgas,j-1),neqgas)
                do m = 1,neqgas
                   do n = 1,neqgas
                      Binv_save(m,n,i,k,j) = ay(j)*Binv(m,n)
                      Bb(m,n,j) = Bb(m,n,j) - &
                           Binv_save(m,n,i,k,j) * cy(j-1)
                   enddo
                enddo
             endif
!             fyv(1:neqgas,j) = fyv(1:neqgas,j) -  &
!                  MATMUL(Binv_save(1:neqgas,1:neqgas,i,k,j),fyv(1:neqgas,j-1))

             do m = 1,neqgas
                ftmp(m) = 0.0d0
                do n =1,neqgas
                   ftmp(m)   = ftmp(m)  + Binv_save(m,n,i,k,j) * fyv(n,j-1)
                enddo
             enddo
             do m = 1,neqgas
                fyv(m,j) = fyv(m,j) - ftmp(m)
             enddo

          enddo

!--- assume zero gradient BC at north to simplify
          j=jmax
          if(.NOT. saved_LU_vars) then
             Bb(1:neqgas,1:neqgas,jmax) = Bb(1:neqgas,1:neqgas,jmax-1) + &
                  Bb(1:neqgas,1:neqgas,jmax)*cy(jmax-1)  
          endif             
          fyv(1:neqgas,jmax) = fyv(1:neqgas,jmax-1) 

!--- START BACKWARD SUBSTITUTION
          if(.NOT. saved_LU_vars) then
!             CALL LUSOLVE(Bb(1:neqgas,1:neqgas,jmax),fyv(1:neqgas,jmax),neqgas)
             CALL LUINV(Bb(1:neqgas,1:neqgas,jmax),neqgas)
             L_save(1:neqgas,1:neqgas,i,k,jmax) = Binv(1:neqgas,1:neqgas)
          endif
          do m = 1,neqgas
             px(i,k,j,m) = 0.0d0
             do n =1,neqgas
                px(i,k,j,m) =  px(i,k,j,m) +  L_save(m,n,i,k,j) * fyv(n,j)
             enddo
          enddo


          do j=jmax-1,js,-1
             OLD_TEMP = px(i,k,j,1)
             fyv(1:neqgas,j) = fyv(1:neqgas,j) - cy(j)* px(i,k,j+1,1:neqgas)
             if(.NOT. saved_LU_vars) then
                CALL LUINV(Bb(1:neqgas,1:neqgas,j),neqgas)
                L_save(1:neqgas,1:neqgas,i,k,j) = Binv(1:neqgas,1:neqgas)
             endif
!             fyv(1:neqgas,j) = MATMUL( L_save(1:neqgas,1:neqgas,i,k,j), fyv(1:neqgas,j) )
             do m = 1,neqgas
                px(i,k,j,m) = 0.0d0
                do n =1,neqgas
                   px(i,k,j,m) =  px(i,k,j,m) +  L_save(m,n,i,k,j) * fyv(n,j)
                enddo
             enddo
!            this evals the error on the last update only
             if(icount == nu)  then 
                wrk_vec(1,i,k,j) = px(i,k,j,1) - OLD_TEMP
             endif
                
          enddo

       ENDDO


          eqn =neqmax
          myres = 0.0d0
          do 15 indx = 1 , eqn_nclr(icolor,0)

             i = eqn_iclr(indx,icolor,0)
             k = eqn_kclr(indx,icolor,0)

             Ts =  px(i,k,0,neqmax)
             fy(js) = Ts
             ay(js) = 0.0d0
             cy(js) = 0.0d0
             by(js) = 1.0d0

             do j = js+1,jmax-1
                ay(j) = - csld(4,i,k,j)
                cy(j) = - csld(5,i,k,j)
                by(j) = csld(1,i,k,j)
                c2 = csld(2,i,k,j)
                c3 = csld(3,i,k,j)
                c8 = csld(8,i,k,j)
                gxy = (+lambdas(i+1,k,j) * (px(i+1,k,j+1,eqn)&
                     - px(i+1,k,j-1,eqn))&
                     - lambdas(i-1,k,j) * (px(i-1,k,j+1,eqn) &
                     - px(i-1,k,j-1,eqn))                    &
                     + lambdas(i,k,j+1) * (px(i+1,k,j+1,eqn) &
                     - px(i-1,k,j+1,eqn))                    &
                     - lambdas(i,k,j-1) * (px(i+1,k,j-1,eqn) &
                     - px(i-1,k,j-1,eqn)))
                fy(j) = dfdt(i,k,j,eqn) + c2*px(i-1,k,j,eqn)  &
                     + c3*px(i+1,k,j,eqn) + c8*gxy
             enddo
             fy(jmax) = tcold
             ay(jmax) = 0.0d0
             cy(jmax) = 0.0d0
             by(jmax) = 1.0d0

! --- THOMAS algorithm
             do j=js+1,jmax
                ay(j) = ay(j)/by(j-1)
                by(j) = by(j)-ay(j)*cy(j-1)
                fy(j) = fy(j)-ay(j)*fy(j-1)
             enddo

             fy(jmax) = fy(jmax)/by(jmax)
             px(i,k,jmax,eqn) = fy(jmax)

             do j=jmax-1,js,-1
                fy(j) = (fy(j)-cy(j)*fy(j+1))/by(j)
!                if(icount == nu)  wrk_vec(2,i,k,j) = px(i,k,j,eqn) - fy(j)
                px(i,k,j,eqn) = fy(j)
             enddo

15       continue     !indx

!          CALL BC_IMP(icolor)  !this calls the interface BC for 2nd order time accuracy

24  continue   !icolor
          saved_LU_vars = .TRUE.
25  continue   !icount
!

          
    col(1:2) = 2
    call parallel_swap4(col,current_mesh)

    IF(resid_flag) then
      myres = zero

!---CHECK convergence
      mysum = zero
      do i=0,xe
      do k=0,ze
      do j=ys,ye        
         mysum(1) = mysum(1) + (wrk_vec(1,i,k,j)/px(i,k,j,1))**2     
         mysum(2) = mysum(2) + (wrk_vec(2,i,k,j)/px(i,k,j,neqmax))**2
      enddo
      enddo
      enddo

!......................EVAL CONVERGENCE.......................................
      if(current_mesh%mesh_num == 1) then
         CALL MPI_ALLREDUCE(mysum,sum,2,MPI_DOUBLE_PRECISION,MPI_SUM,tmpcomm,ierr)
         n_samples = nx*ny
         sum = sqrt(sum/n_samples)
         converged = maxval(sum) < prec_gauss
         maxres(1:2) = sum(1:2)
!         if(myid == 0) write(*,*)icycle,'TXY_LINE',&
!              maxres,converged,current_mesh%mesh_num,prec_gauss
!          write(89,*) maxval(maxres),times(:,1)
      endif
!...........................................................................

    ENDIF   !if (resid_flag)

!---------------------------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE gauss_relax4
!***************************************************************************************************


!***************************************************************************************************
SUBROUTINE gauss_relax5 (current_mesh, nu, resid_flag)

  USE data_types
  USE GLOBAL_DATA 
  USE BC_DATA
  USE LUSOLVER, ONLY : LUSOLVE, LUINV
  USE MYMPI

  IMPLICIT NONE

  include "parallel.h"

!---------------------------------------------------------------------------------------------------
!   Dummy variables:

    TYPE(mesh_pointers) :: current_mesh

    INTEGER, INTENT(IN) :: nu                ! number of relaxations 
    LOGICAL, INTENT(IN) :: resid_flag        ! if this is true it evaluates the residual
   
!   Local Variables
    REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: px
    REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: rhs

    REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: cm,cs

    INTEGER :: icount,i,j,k,xe,ze,ys,ye,m,col(7)
    INTEGER :: indx,icolor,rcolor, tmpcomm
    INTEGER :: eqn1,neqn,eqn,mn,jp1,iteqn,js,jmax, n_samples
    REAL (KIND=double) :: coe,c2,c3,c4,c5,c8
    REAL (KIND=double) :: f1,f2,f3,f4,f5,gxy
    REAL (KIND=double) :: deltaf,addrate,rr,rhsg,myres(2)
    REAL (KIND=double) :: c11,sumc,fx,fz,aa0,a0,sum(2),mysum(2)
!---------------------------------------------------------------------------------------------------

!
!   RETURN immediatly if the processor has no points
!
    if(current_mesh%blank) RETURN

    px => current_mesh%w ; rhs => current_mesh%r
    cm => current_mesh%cm
    cs => csld

    xe = current_mesh%xe ; ze = current_mesh%ze
    ys = 1; ye = ny -1
    eqn1 = 1; neqn=neqmax
    col(5)=0;col(6)=ye;
    col(7) = 1;
    mn = current_mesh%mesh_num

    tmpcomm = current_mesh%comm3d
!--- Start and end of line relaxation loops
    jmax = ny
    js = 0 
!
!   START ITERATIONS
! 
    col(3)=eqn1;col(4)=neqn;
    do 25  icount = 1,nu

    do 25 icolor = 1,2

       rcolor=mod(icolor+2,2)+1     !2,1
       col(1:2) = rcolor
       call parallel_swap4(col,current_mesh)

       do 10 indx = 1 , eqn_nclr(icolor,0)

          i = eqn_iclr(indx,icolor,0)
          k = eqn_kclr(indx,icolor,0)

          Bb(1:neqgas,1:neqgas,js:jmax) = Bm(1:neqgas,1:neqgas,i,k,js:jmax)

          fy(ny) = tcold
          ay(ny) = 0.0d0
          cy(ny) = 0.0d0
          by(ny) = 1.0d0
          zy(ny) = 0.0d0
          eqn = neqmax
          do j = ny-1,1,-1
             ay(j) = - cs(5,i,k,j)
             cy(j) = - cs(4,i,k,j)
             by(j) = cs(1,i,k,j)
             zy(j) = -cs(10,i,k,j)
             c2 = cs(2,i,k,j)
             c3 = cs(3,i,k,j)
             c8 = cs(8,i,k,j)
             gxy = (+lambdas(i+1,k,j) * (px(i+1,k,j+1,eqn)&
                  - px(i+1,k,j-1,eqn))&
                  - lambdas(i-1,k,j) * (px(i-1,k,j+1,eqn) &
                  - px(i-1,k,j-1,eqn))                    &
                  + lambdas(i,k,j+1) * (px(i+1,k,j+1,eqn) &
                  - px(i-1,k,j+1,eqn))                    &
                  - lambdas(i,k,j-1) * (px(i+1,k,j-1,eqn) &
                  - px(i-1,k,j-1,eqn)))

             fy(j) = rhs(i,k,j,eqn) + c2*px(i-1,k,j,eqn)  &
                  + c3*px(i+1,k,j,eqn) + c8*gxy
          enddo

          j=0
          do eqn =1,neqgas
             m = min(eqn,2)
             fx = px(i+1,k,0,eqn)-px(i-1,k,0,eqn)
             fz = px(i,k+1,0,eqn)-px(i,k-1,0,eqn)
             fyv(eqn,j) = capr(eqn,i,k) + capa(m,i,k)*fx+capb(m,i,k)*fz
             Vdiag0(eqn) = capc(m,i,k)
          enddo

! --- Forward sweep solid
          do j=ny-1,1,-1
             ay(j) = ay(j)/by(j+1)
             by(j) = by(j)-ay(j)*cy(j+1)
             fy(j) = fy(j)-ay(j)*fy(j+1)
             zy(j) = zy(j)-ay(j)*zy(j+1)
          enddo
          cy(1) = cy(1) + zy(1)
          zy(1) = 0.0d0

          aa0 = capee(i,k)/by(2)
          a0 = cape(i,k) - aa0*cy(2)
          fyv(1,0) = fyv(1,0) - aa0*fy(2)
          Vdiag0(1) = Vdiag0(1) - aa0*zy(2)
          aa0 = a0/by(1)
          Vdiag0(1) = Vdiag0(1) - aa0*cy(1)
          fyv(1,0) = fyv(1,0) - aa0*fy(1)

!---Fill gas line-coefficients          
          cq(js) = capd(i,k)
          aq(js) = capdd(i,k)     !last term
          do j=ys,ye
             jp1=merge(j+1,j,j<ye)
             c2=cm(2,i,k,j)
             c3=cm(3,i,k,j)
             c8=cm(8,i,k,j)
             aq(j) = - cm(4,i,k,j)
             cq(j) = - cm(5,i,k,j)

             do eqn = 1,neqgas
                f2=px(i-1,k,j,eqn)
                f3=px(i+1,k,j,eqn)
                gxy = (px(i+1,k,jp1,eqn)-px(i+1,k,j-1,eqn)&
                     -px(i-1,k,jp1,eqn)+px(i-1,k,j-1,eqn))

                fyv(eqn,j) = rhs(i,k,j,eqn) + c2*f2+c3*f3 + c8*gxy
             enddo !eqns

          enddo  !j
          fyv(1:neqgas,jmax) = 0.0d0
          aq(jmax) = -1.0d0
          cq(jmax) = 0.0d0

          j = 1
          do eqn = 1,neqgas
             aa0 = aq(j) / Vdiag0(eqn)
             Bb(eqn,eqn,j) = Bb(eqn,eqn,j) - aa0*cq(j-1)
             Vright1(eqn) = cq(j) - aa0*aq(j-1) 
             fyv(eqn,j) = fyv(eqn,j) -  aa0*fyv(eqn,j-1)
          enddo
          
          j = 2
          CALL LUINV(Bb(1:neqgas,1:neqgas,j-1),neqgas)
          Binv = aq(j)*Binv
          do eqn = 1,neqgas
             Bb(1:neqgas,eqn,j) = Bb(1:neqgas,eqn,j) -  Binv(1:neqgas,eqn)*Vright1(eqn)
          enddo
          fyv(1:neqgas,j) = fyv(1:neqgas,j) -  MATMUL(Binv,fyv(1:neqgas,j-1))

          

! --- Forward sweep, GAS
          do j =3,jmax-1
             CALL LUINV(Bb(1:neqgas,1:neqgas,j-1),neqgas)
             Binv = aq(j)*Binv
             Bb(1:neqgas,1:neqgas,j) = Bb(1:neqgas,1:neqgas,j) - Binv*cq(j-1)
             fyv(1:neqgas,j) = fyv(1:neqgas,j) -  MATMUL(Binv,fyv(1:neqgas,j-1))
          enddo

!--- assume zero gradient BC at north to simplify
          Bb(1:neqgas,1:neqgas,jmax) = Bb(1:neqgas,1:neqgas,jmax-1) + &
               Bb(1:neqgas,1:neqgas,jmax)*cq(jmax-1)   
          fyv(1:neqgas,jmax) = fyv(1:neqgas,jmax-1) 

          CALL LUSOLVE(Bb(1:neqgas,1:neqgas,jmax),fyv(1:neqgas,jmax),neqgas)
          px(i,k,jmax,1:neqgas) = fyv(1:neqgas,jmax)

          do j=jmax-1,2,-1
             fyv(1:neqgas,j) = fyv(1:neqgas,j) - cq(j)*fyv(1:neqgas,j+1)
             CALL LUSOLVE(Bb(1:neqgas,1:neqgas,j),fyv(1:neqgas,j),neqgas)
             if(icount == nu)  wrk_vec(1,i,k,j) = px(i,k,j,1) - fyv(1,j)
             px(i,k,j,1:neqgas) = fyv(1:neqgas,j)
          enddo

          j=1
          do eqn=1,neqgas
             fyv(eqn,j) = fyv(eqn,j) - Vright1(eqn)*fyv(eqn,j+1)
          enddo
          CALL LUSOLVE(Bb(1:neqgas,1:neqgas,j),fyv(1:neqgas,j),neqgas)
          px(i,k,j,1:neqgas) = fyv(1:neqgas,j)
 
          j=0          
          do eqn=1,neqgas
             fyv(eqn,j) = (fyv(eqn,j) - cq(j)*fyv(eqn,j+1) - aq(j)*fyv(eqn,j+2))/Vdiag0(eqn)
          enddo
          px(i,k,j,1:neqgas) = fyv(1:neqgas,j)

!--- NOW GO DOWN IN THE SOLID
          fy(0) = px(i,k,0,1)
          px(i,k,0,neqmax) = fy(0)
          do j=1,ny
             fy(j) = (fy(j)-cy(j)*fy(j-1)-zy(j)*fy(0))/by(j)
             if(icount == nu)  wrk_vec(2,i,k,j) = px(i,k,j,neqmax) - fy(j)
             px(i,k,j,neqmax) = fy(j)
          enddo

10     continue   !indx

 
25  continue   !icount,icolor
!
    col(1:2) = 2
    call parallel_swap4(col,current_mesh)
    
    IF(resid_flag) then
       myres = zero

!---CHECK convergence
       mysum = zero
       do i=0,xe
          do k=0,ze
             do j=ys,ye        
                mysum(1) = mysum(1) + (wrk_vec(1,i,k,j)/px(i,k,j,1))**2     
                mysum(2) = mysum(2) + (wrk_vec(2,i,k,j)/px(i,k,j,neqmax))**2
             enddo
          enddo
       enddo


!......................EVAL CONVERGENCE.......................................
       if(current_mesh%mesh_num == 1) then
          CALL MPI_ALLREDUCE(mysum,sum,2,MPI_DOUBLE_PRECISION,MPI_SUM,tmpcomm,ierr)
          n_samples = nx*ny
          sum = sqrt(sum/n_samples)
          converged = maxval(sum) < prec_gauss
          maxres(1:2) = sum(1:2)
!          if(myid == 0) write(*,*)icycle,'TXY_RES',&
!               maxres,converged,current_mesh%mesh_num,myid
       endif
!...........................................................................

    ENDIF   !if (resid_flag)

!---------------------------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE gauss_relax5
!***************************************************************************************************


!***************************************************************************************************
  SUBROUTINE point_relax (current_mesh, nu, resid_flag)

    USE data_types
    USE GLOBAL_DATA 
    USE LUSOLVER, ONLY : LUSOLVE, LUINV
    USE MYMPI

    IMPLICIT NONE

    include "parallel.h"

!---------------------------------------------------------------------------------------------------
!   Dummy variables:

    TYPE(mesh_pointers) :: current_mesh

    INTEGER, INTENT(IN) :: nu                ! number of relaxations 
    LOGICAL, INTENT(IN) :: resid_flag        ! if this is true it evaluates the residual

!   Local Variables
    REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: px
    REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: rhs

    REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: cm,cs

    INTEGER :: icount,i,j,k,xe,ze,ys,ye,m,col(7)
    INTEGER :: indx,icolor,rcolor, tmpcomm
    INTEGER :: eqn1,neqn,eqn,mn,jp1,iteqn,js,jmax,n_samples
    REAL (KIND=double) :: coe,c2,c3,c4,c5,c8,Ts
    REAL (KIND=double) :: f1,f2,f3,f4,f5,gxy
    REAL (KIND=double) :: deltaf,addrate,rr,rhsg,myres(2),resmax
    REAL (KIND=double) :: c11,sumc,omg_txy,one_m_omg, sum(2), mysum(2)
!---------------------------------------------------------------------------------------------------

!
!   RETURN immediatly if the processor has no points
!
    if(current_mesh%blank) RETURN

    px => current_mesh%w ; rhs => current_mesh%r
    cm => current_mesh%cm
    cs => csld

    xe = current_mesh%xe ; ze = current_mesh%ze
    ys = 1; ye = ny -1
    eqn1 = 1; neqn=neqmax
    col(5)=0;col(6)=ye;
    col(7) = 1
    mn = current_mesh%mesh_num

    tmpcomm = current_mesh%comm3d
    omg_txy=1.00d0
    one_m_omg = (one - omg_txy)

!--- Start and end of line relaxation loops
    jmax = ny
    js = 0 
!
!   START ITERATIONS
! 
    col(3)=eqn1;col(4)=neqn;
    do   icount = 1,nu


       do  icolor = 1,2


          rcolor=mod(icolor+2,2)+1     !2,1
          col(1:2) = rcolor
          call parallel_swap4(col,current_mesh)



          do  indx = 1 , eqn_nclr(icolor,0)

             i = eqn_iclr(indx,icolor,0)
             k = eqn_kclr(indx,icolor,0)

             do j=ys,ye
                if(j<ye)then
                   jp1 = j+1
                else
                   jp1 = j
                endif
                c2 = cm(2,i,k,j)
                c3 = cm(3,i,k,j)
                c4 = cm(4,i,k,j)
                c5 = cm(5,i,k,j)
                c8 = cm(8,i,k,j)

                do eqn = 1,neqgas
                   f2=px(i-1,k,j,eqn)
                   f3=px(i+1,k,j,eqn)
                   f4=px(i,k,j-1,eqn)
                   f5=px(i,k,j+1,eqn)

                   gxy = (px(i+1,k,jp1,eqn)-px(i+1,k,j-1,eqn)&
                        -px(i-1,k,jp1,eqn)+px(i-1,k,j-1,eqn))

                   addrate=zero
                   do m=1,4
                      if(m.ne.eqn) addrate = addrate - Bm(eqn,m,i,k,j)*px(i,k,j,m)
                   enddo
                   rr = rhs(i,k,j,eqn)+addrate
                   rhsg = c2*f2+c3*f3+c4*f4+c5*f5+c8*gxy+rr
                   coe = 1.0d0/Bm(eqn,eqn,i,k,j)
                   deltaf = rhsg*coe - px(i,k,j,eqn)

                   if(icount == nu)  wrk_vec(1,i,k,j) = deltaf
                   px(i,k,j,eqn) =  px(i,k,j,eqn) + omg_txy*deltaf

                enddo   !eqns

             enddo  !j 


          enddo  !indx


          eqn =neqmax
          myres = 0.0d0
          do  indx = 1 , eqn_nclr(icolor,0)

             i = eqn_iclr(indx,icolor,0)
             k = eqn_kclr(indx,icolor,0)

             Ts =  px(i,k,0,neqmax)
             fy(js) = Ts
             ay(js) = 0.0d0
             cy(js) = 0.0d0
             by(js) = 1.0d0
             do j = js+1,jmax-1
                ay(j) = - cs(4,i,k,j)
                cy(j) = - cs(5,i,k,j)
                by(j) = cs(1,i,k,j)
                c2 = cs(2,i,k,j)
                c3 = cs(3,i,k,j)
                c8 = cs(8,i,k,j)
                gxy = (+lambdas(i+1,k,j) * (px(i+1,k,j+1,eqn)&
                     - px(i+1,k,j-1,eqn))&
                     - lambdas(i-1,k,j) * (px(i-1,k,j+1,eqn) &
                     - px(i-1,k,j-1,eqn))                    &
                     + lambdas(i,k,j+1) * (px(i+1,k,j+1,eqn) &
                     - px(i-1,k,j+1,eqn))                    &
                     - lambdas(i,k,j-1) * (px(i+1,k,j-1,eqn) &
                     - px(i-1,k,j-1,eqn)))
                fy(j) = rhs(i,k,j,eqn) + c2*px(i-1,k,j,eqn)  &
                     + c3*px(i+1,k,j,eqn) + c8*gxy
             enddo
             fy(jmax) = tcold
             ay(jmax) = 0.0d0
             cy(jmax) = 0.0d0
             by(jmax) = 1.0d0

! --- THOMAS algorithm
             do j=js+1,jmax
                ay(j) = ay(j)/by(j-1)
                by(j) = by(j)-ay(j)*cy(j-1)
                fy(j) = fy(j)-ay(j)*fy(j-1)
             enddo

             fy(jmax) = fy(jmax)/by(jmax)
             px(i,k,jmax,eqn) = fy(jmax)

             do j=jmax-1,js,-1
                fy(j) = (fy(j)-cy(j)*fy(j+1))/by(j)
                if(icount == nu)  wrk_vec(2,i,k,j) = px(i,k,j,eqn) - fy(j)
                px(i,k,j,eqn) = fy(j)
             enddo

          enddo

!          CALL BC_IMP(icolor)

       enddo   !icolor
    enddo   !icount
!
    col(1:2) = 2
    call parallel_swap4(col,current_mesh)

    IF(resid_flag) then
       myres = zero

!---CHECK convergence
       mysum = zero
       do i=0,xe
          do k=0,ze
             do j=ys,ye        
                mysum(1) = mysum(1) + (wrk_vec(1,i,k,j)/px(i,k,j,1))**2     
                mysum(2) = mysum(2) + (wrk_vec(2,i,k,j)/px(i,k,j,neqmax))**2
             enddo
          enddo
       enddo

!......................EVAL CONVERGENCE.......................................
       if(current_mesh%mesh_num == 1) then
          CALL MPI_ALLREDUCE(mysum,sum,2,MPI_DOUBLE_PRECISION,MPI_SUM,tmpcomm,ierr)
          n_samples = nx*ny
          sum = sqrt(sum/n_samples)
!          converged = maxval(sum) < prec_gauss
          maxres(1:2) = sum(1:2)
         if(myid == 0) write(*,*)icycle,'POINT_RELAX_TXY',&
              maxres,converged,current_mesh%mesh_num,prec_gauss
          write(39,*) maxval(maxres),times(:,1)
       endif

!...........................................................................

    ENDIF   !if (resid_flag)
!---------------------------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE point_relax
!***************************************************************************************************
