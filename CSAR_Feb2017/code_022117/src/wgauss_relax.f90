MODULE IMPLICIT
  Use Timing

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
!
  REAL*8,PARAMETER :: SMALL = 1d-30
  Logical :: implicitBC =.TRUE.
  INTEGER :: iterLINEARsys ,iterativemethod ,GAUSS_GMRES,col(7)
  INTEGER :: GMRES_restarts,NewtonIteration
  REAL*8 ::  TOLERANCEGAUSS = 1d-12


!***************************************************************************************************

contains

!****************************************************************************************
  SUBROUTINE ADJUST_LHS

    USE GLOBAL_DATA

    IMPLICIT NONE

    include 'parallel.h'

!----------------------------------------------------------------------------------------
! Local variables
    TYPE(mesh_pointers), POINTER :: current_mesh
    INTEGER ::  i,j,k,ii,jj,kk, mn, ys,ye,xs,xe,zs,ze
    INTEGER ::  m,n, ip,kp,jp, col(4)
    REAL*8  :: corrflux,s1,s2,s3,sig1,sig2,dx1a,dz1a,dy1a,dya,minmod
    REAL*8  :: nu,nux,nuz,nuy
!----------------------------------------------------------------------------------------

!
    current_mesh => finest_mesh
    col(1) = 1; col(2) = neqgas; col(3) = 0; col(4) = drange(6)
    if (it_new > 1 ) CALL TWOLAYER_SWAP(col,current_mesh)

    ys = 1
    ye = drange(6)-1
    xs = drange(1)
    xe = drange(2)
    zs = drange(3)
    ze = drange(4)

!
! N-S Boundary cond for B
!
    Bm(1:neqgas,1:neqgas,:,:,0) = 0.0d0
    Bm(1:neqgas,1:neqgas,:,:,ny) = 0.0d0
    do m = 1,neqgas
       Bm(m,m,:,:,0) = 1.0d0
       Bm(m,m,:,:,ny) = 1.0d0
    enddo

    DO j = ys, ye
       do k = zs, ze
          do i = xs, xe
             do m=1,neqgas
                Bm(m,m,i,k,j) = coeff(m,1,i,k,j) + Bm(m,m,i,k,j)          !for line relaxation
             enddo
          enddo
       enddo
    ENDDO

!------------------------------------------------------------------------------------
    RETURN 
  END SUBROUTINE ADJUST_LHS
!****************************************************************************************

  SUBROUTINE gauss_relax4 (current_mesh, r4, nu, resid_flag)

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
    real*8,  Intent(INOUT) :: r4(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize)
    LOGICAL, INTENT(IN) :: resid_flag        ! if this is true it evaluates the residual

!   Local Variables
!!>    REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: px
!!>    REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: rhs


    INTEGER :: icount,i,j,k,xe,ze,ys,ye,m,n,jj,jjj
    INTEGER :: indx,icolor,rcolor, tmpcomm,wunit,idir
    INTEGER :: eqn1,neqn,eqn,mn,jp1,iteqn,js,jmax,ierrWG
    INTEGER,DIMENSION(3) :: ikj,ikj_P, ikj_P1,ikj_M1,ikj_M2,ikjP,ikjM
    REAL (KIND=double) :: coe,c2,c3,c6,c7,c4,c5,c8,c9,Ts
    REAL (KIND=double) :: f1,f2,f3,f6,f7,f4,f5,gxy,gzy
    REAL (KIND=double) :: deltaf,addrate,rr,rhsg,myres(2),resmax
    REAL (KIND=double) :: c11,sumc,omg_txy, sum(2), mysum(2)
    REAL (KIND=double) :: ftmp(maxsize), OLD_TEMP, corrflux(maxsize),nuy
    REAL (KIND=double) :: s1,s2,s3,sig1,sig2
    LOGICAL outcome,repeat
!---------------------------------------------------------------------------------------------------


    wunit =6  !for writing errors
    qpick5 = (/1,0,0, 0, 5/)
!
!   RETURN immediatly if the processor has no points
!
    if(current_mesh%blank) RETURN

!!>    px => current_mesh%w ; rhs => current_mesh%r

    xe = current_mesh%xe ; ze = current_mesh%ze
    ys = 1; ye = ny -1
    eqn1 = 1; neqn=maxsize
    col(5)=0;col(6)=ye+1;
    col(7) = 1
    mn = current_mesh%mesh_num

    tmpcomm = current_mesh%comm3d
    omg_txy=1.00

!--- Start and end of line relaxation loops

    js = 0
!BC
    fyv(1:neqgas,:) = zero
    fy(:) = zero
    ay(js) = 0d0
    by(js) = 1d0
    cy(js) = 0d0

!
!   START ITERATIONS
! 
    col(3)=eqn1;col(4)=neqn;
    gassolidloop:do icount = 1,nu

       colorsloop:do icolor = 1,2


!fluid part
          if(JRFLO(1) >=ny) then
             jmax = ye+1
             ay(jmax) = -1.0d0
          else
             jmax = JRFLO(1)
             ay(jmax) = 0.0d0
          endif
          cy(jmax) = 0.0d0
          fyv(1:neqgas,jmax) = zero

          rcolor=mod(icolor+2,2)+1     !2,1
          col(1:2) = rcolor
          call parallel_swap4(col,current_mesh)

          if(implicitBC) then
             CALL BC_IMP(icolor,r4,Vfun = dff)
          endif
          call gauss_solid
          if(implicitBC) then
             CALL BC_IMP(icolor,r4,Vfun = dff)
          endif
          call gauss_gas  
          if(implicitBC) then
             CALL BC_IMP(icolor,r4,Vfun = dff)
          endif
          call gauss_solid


       enddo colorsloop
       saved_LU_vars = .TRUE.
    enddo gassolidloop
!

    col(1:2) = 2
    call parallel_swap4(col,current_mesh)

    IF(resid_flag) then
       myres = zero

       eqn =maxsize
       do i = drange(1),drange(2)
          do k = drange(3),drange(4)
             do j = drange(5),drange(6)
                c2=csld(2,i,k,j)
                c3=csld(3,i,k,j)
                c6=csld(6,i,k,j)
                c7=csld(7,i,k,j)
                c4=csld(4,i,k,j)
                c5=csld(5,i,k,j)
                c8=csld(8,i,k,j)
                c9=csld(9,i,k,j)
                coe = csld(1,i,k,j)


                f2=dff(i-1,k,j,eqn)
                f3=dff(i+1,k,j,eqn)
                f6=dff(i,k,j-1,eqn)
                f7=dff(i,k,j+1,eqn)
                f4=dff(i,k-1,j,eqn)
                f5=dff(i,k+1,j,eqn)
                gxy = (+lambdas(i+1,k,j) * (dff(i+1,k,j+1,eqn)&
                     - dff(i+1,k,j-1,eqn))&
                     -lambdas(i-1,k,j) * (dff(i-1,k,j+1,eqn)&
                     - dff(i-1,k,j-1,eqn))&
                     +lambdas(i,k,j+1) * (dff(i+1,k,j+1,eqn)&
                     - dff(i-1,k,j+1,eqn))&
                     -lambdas(i,k,j-1) * (dff(i+1,k,j-1,eqn)&
                     - dff(i-1,k,j-1,eqn)))
                gzy = (+lambdas(i,k+1,j) * (dff(i,k+1,j+1,eqn)&
                     - dff(i,k+1,j-1,eqn))&
                     -lambdas(i,k-1,j) * (dff(i,k-1,j+1,eqn)&
                     - dff(i,k-1,j-1,eqn))&
                     +lambdas(i,k,j+1) * (dff(i,k+1,j+1,eqn)&
                     - dff(i,k-1,j+1,eqn))&
                     -lambdas(i,k,j-1) * (dff(i,k+1,j-1,eqn)&
                     - dff(i,k-1,j-1,eqn)))
                rr = r4(i,k,j,eqn)
                rhsg = rr+c2*f2+c3*f3+c4*f4+c5*f5+c6*f6+c7*f7+c8*gxy+c9*gzy &
                     - coe*dff(i,k,j,eqn) 

             enddo
          enddo
       enddo
       if(ncyc_run < 100 .OR. n_samples == 0) &   !just pick the first 100 cycleas to be sure
            CALL MPI_ALLREDUCE(m,n_samples,1,MPI_INTEGER,MPI_SUM,tmpcomm,ierrWG)

!......................EVAL CONVERGENCE.......................................
       if(current_mesh%mesh_num == 1) then
          CALL MPI_ALLREDUCE(mysum,sum,2,MPI_DOUBLE_PRECISION,MPI_SUM,tmpcomm,ierr)
          sum = sqrt(sum/n_samples)
          converged = maxval(sum) < prec_gauss
          maxres(1:2) = sum(1:2)
       endif

    ENDIF   !if (resid_flag)
!---------------------------------------------------------------------------------------------------
    RETURN

  contains

    subroutine gauss_Gas

      integer :: indx,i,k,j,l,m,eqn
      real*8 :: omg_txy,one_m_omg 
!---------------------------------------

      jmax = ny
      js = 0
!dirichlet BC
      aq(:,js) = 0.0d0
      cq(:,js) = 0.0d0
      fyv(1:neqgas,js) = 0.0d0

      aq(:,jmax) = -1.0d0
      cq(:,jmax) = 0.0d0
      gasloop:do indx = 1, current_mesh%nclr(icolor)
         omg_txy = 1.00d0
         one_m_omg = one - omg_txy

         i = current_mesh%iclr(indx,icolor)
         k = current_mesh%kclr(indx,icolor)

         if(.NOT. saved_LU_vars) then
            do j = 0,jmax
               do m = 1,neqgas
                  Bb(m,1:neqgas,j) = Bm(m,1:neqgas,i,k,j)
               enddo
            enddo
         endif

         ikj = (/i,k,0/)   !j is set later
         do j=js+1,jmax-1

            if(j<jmax-1)then
               jp1 = j+1
            else
               jp1 = j
            endif
            ikj(3) = j
            do eqn = 1,neqgas
               c2=coeff(eqn,2,i,k,j)
               c3=coeff(eqn,3,i,k,j)
               c4=coeff(eqn,4,i,k,j)
               c5=coeff(eqn,5,i,k,j)
               c8=coeff(eqn,8,i,k,j)
               c9=coeff(eqn,9,i,k,j)
               aq(eqn,j) = - coeff(eqn,6,i,k,j)
               cq(eqn,j) = - coeff(eqn,7,i,k,j)

               f2=dff(i-1,k,j,eqn)
               f3=dff(i+1,k,j,eqn)
               f4=dff(i,k-1,j,eqn)
               f5=dff(i,k+1,j,eqn)
               gxy = (dff(i+1,k,jp1,eqn)-dff(i+1,k,j-1,eqn)&
                    -dff(i-1,k,jp1,eqn)+dff(i-1,k,j-1,eqn))
               gzy = (dff(i,k+1,jp1,eqn)-dff(i,k+1,j-1,eqn)&
                    -dff(i,k-1,jp1,eqn)+dff(i,k-1,j-1,eqn))
               corrflux(eqn) = c2*f2+c3*f3 + c4*f4+c5*f5 &
                    + c8*gxy+c9*gzy 
            enddo
!!>            if(j==1) corrflux(1:neqgas) = corrflux(1:neqgas) - aq(1:neqgas,j)*dff(i,k,j-1,1:neqgas)


            fyv(1:neqgas,j) = r4(i,k,j,1:neqgas) + corrflux(1:neqgas)

         enddo  !j


! --- THOMAS algorithm--start
         do j=js+1,jmax-1
            if(.NOT. saved_LU_vars) then
               CALL LUINV(Bb(1,1,j-1),ubound(Bb,1),neqgas)
               do m = 1,neqgas
                  do n = 1,neqgas
                     Binv_save(m,n,i,k,j) = aq(m,j)*Binv(m,n)
                     Bb(m,n,j) = Bb(m,n,j) - Binv_save(m,n,i,k,j)*cq(n,j-1)
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

         j=jmax
!--- assume zero gradient BC at north to simplify
         if(jmax == ny ) then
            if(.NOT. saved_LU_vars) then
               do n = 1,neqgas
                  do m = 1,neqgas
                     Bb(m,n,jmax) = Bb(m,n,jmax-1) + &
                          cq(m,jmax-1)*deltaKron(m,n)
                  end do
               end do
            endif
            fyv(1:neqgas,jmax) = fyv(1:neqgas,jmax-1) 

!--- START BACKWARD SUBSTITUTION
            if(.NOT. saved_LU_vars) then
!----CALL LUSOLVE(Bb(1,1,jmax),fyv(1:neqgas,jmax),ubound(Bb,1),neqgas)
               CALL LUINV(Bb(1,1,jmax),ubound(Bb,1),neqgas)
               L_save(1:neqgas,1:neqgas,i,k,jmax) = Binv(1:neqgas,1:neqgas)
            endif
            do m = 1,neqgas
               dff(i,k,j,m) = 0.0d0
               do n =1,neqgas
                  dff(i,k,j,m) =  dff(i,k,j,m) +  L_save(m,n,i,k,j) * fyv(n,j)
               enddo
            enddo
         else
            do m = 1,neqgas
               dff(i,k,j,m) = zero
            enddo
         endif

!backward substitution
         jloop: do j=jmax-1,js,-1
            fyv(1:neqgas,j) = fyv(1:neqgas,j) - cq(1:neqgas,j)* dff(i,k,j+1,1:neqgas)
            if(.NOT. saved_LU_vars) then
               CALL LUINV(Bb(1,1,j),ubound(Bb,1),neqgas)
               L_save(1:neqgas,1:neqgas,i,k,j) = Binv(1:neqgas,1:neqgas)
            endif
!             fyv(1:neqgas,j) = MATMUL( L_save(1:neqgas,1:neqgas,i,k,j), fyv(1:neqgas,j) )
            do m = 1,neqgas
               dff(i,k,j,m) = 0.0d0
               do n =1,neqgas
                  dff(i,k,j,m) =  dff(i,k,j,m) +  L_save(m,n,i,k,j) * fyv(n,j)
               enddo
            enddo

            if(ncyc_run < ncyc_debug) then
               call NAN_VAL(f(i,k,j,1:neqgas)+dff(i,k,j,1:neqgas),outcome,neqgas)
               call MPI_ALLREDUCE(outcome,repeat,1,MPI_LOGICAL,MPI_LOR,comm3d,ierrWG)
               if (repeat .and. .not. outcome) then
                  call vizPACK3D(.true.)
                  call MPI_Abort(MPI_COMM_WORLD,ierrWG)
                  call MPI_finalize(m)
                  stop
               elseif(outcome) then
!!>               if(outcome) then
                  jj = j
                  write(wunit,*)'DEBUG SECTION FOR TEMPERATURE AND SPECIES SOLVER'
                  write(wunit,*)'N_CYCLE',ncyc
                  write(wunit,*)'Myid',myid
                  write(wunit,*)'ITER',icount,icolor
                  write(wunit,*)'POS REL This proc(ikj)',i,k,jj
                  write(wunit,*)'Burn rate at this position on surf,psi,temp,j=0',rb(i,k),psi(i,k,0),f(i,k,0,1:maxsize)
                  write(wunit,*)'JRFLO',JRFLO
                  write(wunit,*)'Grid limits-local',drange(2),drange(4),ny
                  write(wunit,*)'POS ABS',i+ib_mg(1),k+kb_mg(1),jj
                  write(wunit,'(a,1p50e12.4)')'NEW SOLN',f(i,k,jj,1:neqgas)+dff(i,k,jj,1:neqgas)
                  write(wunit,'(a,1p50e12.4)')'OLD SOLN',oldsoln(i,k,jj,1:neqgas)
                  write(wunit,*)'New/old solution by j'
                  do jjj = 0,ny
                     write(wunit,'(i3,1p50e12.4)') jjj,f(i,k,jjj,1:neqgas)+dff(i,k,jjj,1:neqgas)
                     write(wunit,'(a3,1p50e12.4)') 'OLD',oldsoln(i,k,jjj,1:neqgas)
                  enddo
                  write(wunit,*)'rate by j'
                  do jjj = 0,ny
                     write(wunit,'(i3,1p50e12.4)') jjj,rate(i,k,jjj,1:neqgas)
                  enddo
                  write(wunit,*)'radiation and properties by j'
                  do jjj = 0,ny
                     write(wunit,'(i3,1p50e12.4)') jjj,dtorhocp(i,k,jjj),&
  radheat(i,k,jjj),lambdag(i,k,jjj),rhog(i,k,jjj),f(i,k,jjj,1),f(i,k,jjj,5:6)
                  enddo
                  write(wunit,'(a,1p50e12.4)')'Explicit', r4(i,k,jj,1:neqgas)
                  write(wunit,'(a,1p50e12.4)')'Corrective', dconv(i,k,jj,1:ndim)
                  write(wunit,*)'BM by row'
                  do m =1,neqgas
                     write(wunit,'(i3,1p50e12.4)')m,Bm(m,1:neqgas,i,k,j), Bb(m,1:neqgas,j)
                  enddo
                  do jjj = 0,ny
                     write(wunit,*)'Coeffs by equation',jjj
                     do m =1,neqgas
                        write(wunit,'(i3,1p50e12.4)')m,coeff(m,1:9,i,k,jjj)
                     enddo
                  enddo
                  call vizPACK3D(.true.)
                  call MPI_Abort(MPI_COMM_WORLD,ierrWG)
                  call MPI_finalize(ierr)
                  stop
               endif
            endif
         enddo jloop

      enddo gasloop
      

      return

    end subroutine gauss_Gas
!**********************************************
!**********************************************
    subroutine gauss_solid
      integer :: indx,i,k,j,m,eqn
!---------------------------------------

      eqn =maxsize
      jmax = ny
      ay(jmax) = 0.0d0
      cy(jmax) = 0.0d0
      by(jmax) = 1.0d0
      myres = 0.0d0
      solidloop:do  indx = 1 , current_mesh%nclr(icolor)

         i = current_mesh%iclr(indx,icolor)
         k = current_mesh%kclr(indx,icolor)


         do j = js+1,jmax-1
            ay(j) = - csld(6,i,k,j)
            cy(j) = - csld(7,i,k,j)
            by(j) = csld(1,i,k,j)
            c2 = csld(2,i,k,j)
            c3 = csld(3,i,k,j)
            c4 = csld(4,i,k,j)
            c5 = csld(5,i,k,j)
            c8 = csld(8,i,k,j)
            c9 = csld(9,i,k,j)
            gxy = (+lambdas(i+1,k,j) * (dff(i+1,k,j+1,eqn)&
                 - dff(i+1,k,j-1,eqn))&
                 - lambdas(i-1,k,j) * (dff(i-1,k,j+1,eqn) &
                 - dff(i-1,k,j-1,eqn))                    &
                 + lambdas(i,k,j+1) * (dff(i+1,k,j+1,eqn) &
                 - dff(i-1,k,j+1,eqn))                    &
                 - lambdas(i,k,j-1) * (dff(i+1,k,j-1,eqn) &
                 - dff(i-1,k,j-1,eqn)))
            gzy = (+lambdas(i,k+1,j) * (dff(i,k+1,j+1,eqn)&
                 - dff(i,k+1,j-1,eqn))&
                 -lambdas(i,k-1,j) * (dff(i,k-1,j+1,eqn)&
                 - dff(i,k-1,j-1,eqn))&
                 +lambdas(i,k,j+1) * (dff(i,k+1,j+1,eqn)&
                 - dff(i,k-1,j+1,eqn))&
                 -lambdas(i,k,j-1) * (dff(i,k+1,j-1,eqn)&
                 - dff(i,k-1,j-1,eqn)))
            fy(j) = r4(i,k,j,eqn) + c2*dff(i-1,k,j,eqn)  &
                 + c3*dff(i+1,k,j,eqn) + c4*dff(i,k-1,j,eqn) + c5*dff(i,k+1,j,eqn) &
                 + c8*gxy + c9*gzy + csld(10,i,k,j)*dff(i,k,0,maxsize)
         enddo

         fy(1) = fy(1) - ay(1)*dff(i,k,0,1)

! --- THOMAS algorithm
         do j=js+1,jmax
            ay(j) = ay(j)/by(j-1)
            by(j) = by(j)-ay(j)*cy(j-1)
            fy(j) = fy(j)-ay(j)*fy(j-1)
         enddo

         fy(jmax) = fy(jmax)/by(jmax)
         dff(i,k,jmax,eqn) = fy(jmax)

         do j=jmax-1,js+1,-1
            fy(j) = (fy(j)-cy(j)*fy(j+1))/by(j)
            dff(i,k,j,eqn) = fy(j)
         enddo

      enddo solidloop
      return
    END SUBROUTINE gauss_solid


  end SUBROUTINE gauss_relax4
!***************************************************************************************************


!***************************************************************************************************
  subroutine GMres_solve(itGMRES,itGS,restarted)
!
!GMRES with SOR as a left preconditioner.
!itGMRES is the max num of GMres iters,itGS is the max numb of preconditioning iters
!
    use global_data
    USE MYMPI

    implicit none
    include 'parallel.h'

    Integer,Intent(IN) :: itGMRES,itGS
    logical,Intent(IN) :: restarted
!local variables
    Integer :: i,k,j,M,it,ii,kk,jj
    real*8 :: betaGMRES,alphaGMRES,htmp,rd,hd,rr,initialresid
    logical :: stoppage
    character(len =5) :: prec
!---------------------------------------

    col(5)=drange(5);col(6)=drange(6);
    col(7) = 1
    col(3)=1;col(4)=maxsize
    col(1) = 1
    col(2) = 2
    prec = 'left '

    dqdt = dff
    restartif: if(restarted) then
       call matvec(dff(-2,-2,-2,1),wrk4)  
       do i = drange(1),drange(2)
          do k = drange(3),drange(4)
             do j = drange(5),drange(6)
                wrk4(i,k,j,:) = dfdt(i,k,j,:) - wrk4(i,k,j,:)
             enddo
          enddo
       enddo
       call vecvec(wrk4,wrk4,betaGMRES)
       betaGMRES = max(sqrt(betaGMRES),SMALL)
       if(prec == 'left ') then
          dff = 0d0
          call gauss_relax4(finest_mesh,wrk4,itGS,.FALSE.) 
       else
          dff = wrk4
       endif
    else
       if(prec == 'left ') then
          call gauss_relax4(finest_mesh,dfdt,itGS,.FALSE.)   !p initial guess should be =0
       else
          dff = dfdt
          call parallel_swap4(col,finest_mesh)
       endif
    endif restartif
!
    call vecvec(dff,dff,betaGMRES)
    betaGMRES = max(sqrt(betaGMRES),SMALL)
    do M = 1,maxsize
       do j = drange(5),drange(6)
          do k = drange(3)-1,drange(4)+1
             do i = drange(1)-1,drange(2)+1
                vgm(i,k,j,M,1) = dff(i,k,j,M)/betaGMRES
             enddo
          enddo
       enddo
!!>       vgm(drange(1)-1:drange(2)+1,drange(3)-1:drange(4)+1,drange(6),M,1) = vgm(drange(1)-1:drange(2)+1,drange(3)-1:drange(4)+1,drange(6)-1,M,1)
    enddo

!
    dgm(1) = betaGMRES;
    initialresid = betaGMRES
    call vecvec(dfdt,dfdt,alphaGMRES)
    alphaGMRES = max(sqrt(alphaGMRES),SMALL)

    kk = 0
    GmresLOOP:do while (.true.)
       kk = kk+1
       stoppage = kk >= min(ubound(vgm,5),itGMRES)-1


       if(Prec == 'right') then
          dff = 0d0
          call gauss_relax4(finest_mesh,vgm(-2,-2,-2,1,kk),itGS,.FALSE.)
          call matvec(dff,wrk4)
       else
          call matvec(vgm(-2,-2,-2,1,kk),wrk4)   !it needs ghost points but not BC points
          if(prec == 'left ') then
             dff = zero
             call gauss_relax4(finest_mesh,wrk4,itGS,.FALSE.) 
          else
             dff = wrk4
             call parallel_swap4(col,finest_mesh)
          endif
       endif

       do ii = 1,kk
          call vecvec(dff,vgm(-2,-2,-2,1,ii),hgm(ii,kk))
       enddo

       do ii = 1,kk
          do m =1,maxsize
             do j = drange(5),drange(6)
                do k = drange(3)-1,drange(4)+1
                   do i = drange(1)-1,drange(2)+1
                      dff(i,k,j,m) = dff(i,k,j,m) - hgm(ii,kk)*vgm(i,k,j,m,ii)
                   enddo
                enddo
             enddo
          enddo
       enddo

       call vecvec(dff,dff,betaGMRES)
       betaGMRES = max(sqrt(betaGMRES),SMALL)
       hgm(kk+1,kk) = betaGMRES

       do M = 1,maxsize
          do j = drange(5),drange(6)
             do k = drange(3)-1,drange(4)+1
                do i = drange(1)-1,drange(2)+1
                   vgm(i,k,j,M,kk+1) = dff(i,k,j,M)/betaGMRES
                enddo
             enddo
          enddo
!!>          vgm(drange(1)-1:drange(2)+1,drange(3)-1:drange(4)+1,drange(6),M,kk+1) = vgm(drange(1)-1:drange(2)+1,drange(3)-1:drange(4)+1,drange(6)-1,M,kk+1)
       enddo

       do ii = 1,kk-1
          htmp = c1gm(ii) * hgm(ii,kk) + c2gm(ii) * hgm(ii+1,kk);
          hgm(ii+1,kk) = - c2gm(ii) * hgm(ii,kk) + c1gm(ii) * hgm(ii+1,kk);
          hgm(ii,kk) = htmp;
       end do
       rD = hgm(kk,kk);
       hD = hgm(kk+1,kk);
       rr = max(sqrt(rD*rD+hD*hD),SMALL)
       c1gm(kk) = rD/rr;
       c2gm(kk) = hD/rr;
       hgm(kk,kk) = rr;
       hgm(kk+1,kk) = zero;
       dgm(kk+1) = -c2gm(kk)*dgm(kk);
       dgm(kk) = c1gm(kk)*dgm(kk);
       errGMres = abs(dgm(kk+1))!/initialresid
!!>       if (myid == 0) write(6,*)kk,'errGMres =',errGMres
       if (errGMres  < TOLERANCEGAUSS .or. stoppage) then
!-----          if(myid == 0)write(*,*)'GMresAbs',abs(dgm(kk+1))
          do ii = kk,1,-1
             ygm(ii) = dgm(ii)/hgm(ii,ii);
             do jj = ii-1,1,-1
                dgm(jj) = dgm(jj) - hgm(jj,ii)*ygm(ii);
             end do
             if(prec == 'right') then
                dff = 0d0
                call gauss_relax4(finest_mesh,vgm(-2,-2,-2,1,ii),itGS,.FALSE.)
                dqdt = dqdt + ygm(ii)*dff
             else
                do M = 1,maxsize
                   do j = drange(5),drange(6)
                      do k = drange(3)-1,drange(4)+1
                         do i = drange(1)-1,drange(2)+1
                            dqdt(i,k,j,M) = dqdt(i,k,j,M) + ygm(ii)*vgm(i,k,j,M,ii) 
                         enddo
                      enddo
                   enddo
                enddo
             endif
          end do
          exit GmresLOOP !done
       endif
    enddo GmresLOOP
    
    dff = dqdt
    numerdivergence = .not. errGMres < 1d99
    

!-----    write(*,*) dff(0,0,0,1),dff(0,0,0,maxsize)
!-----    call matvec(dff(-2,-2,-2,1),wrk4)  
!-----    do i = drange(1),drange(2)
!-----       do k = drange(3),drange(4)
!-----          do j = drange(5),drange(6)
!-----             wrk4(i,k,j,:) = dfdt(i,k,j,:) - wrk4(i,k,j,:)
!-----          enddo
!-----       enddo
!-----    enddo
!-----    if(implicitBC)write(*,'(a,1p5e12.5,a,1p5e12.5)')'Surf REsid',wrk4(0,0,0,1:maxsize)
!-----    do j=0,ny
!-----       write(*,'(i3,1p10e12.4)')j,wrk4(0,0,j,1:maxsize)
!-----    end do
!-----    call vecvec(wrk4,wrk4,betaGMRES)
!-----    betaGMRES = sqrt(betaGMRES)
!-----    if(myid == 0) write(*,*)'Final residual not precond',betaGMRES
!-----    dff = 0d0
!-----    call gauss_relax4(finest_mesh,wrk4,itGS,.FALSE.) 
!-----    do j=0,ny
!-----       write(*,'(i3,1p10e12.4)')j,dff(0,0,j,1:maxsize)
!-----    end do
!-----    call vecvec(dff,dff,betaGMRES)
!-----    betaGMRES = sqrt(betaGMRES)
!-----    if(myid == 0) write(*,*)'Final residual  precond',betaGMRES
!-----    stop

    return
  end subroutine gmres_solve
!----------------------------

!----------------------------
  subroutine matvec(vIN,vOUT)

    use global_data
    USE MYMPI 
    USE BC_DATA

    implicit none

    real*8,Intent(INOUT) ::  vIN (-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize)
    real*8,Intent(INOUT) :: vOUT(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize)
!local variables
    Integer :: i,k,j,M,it,jp1,eqn
    Real*8 :: tmp(maxsize),mynormresid
    Real*8 :: c2,c3,c6,c7,c4,c5,c8,c9,f2,f1,f3,f6,f7,f4,f5,gxy 
    Real*8 :: gzy,addrate,rr,coe,rhsg,fx,fz,term
!---------------------------------------


    do eqn=1,neqgas
       do j = drange(5)+1,drange(6)-1
          do k = drange(3),drange(4)
             do i = drange(1),drange(2)
                jp1=merge(j+1,j,j<drange(6)-1)
                c2=coeff(eqn,2,i,k,j)
                c3=coeff(eqn,3,i,k,j)
                c6=coeff(eqn,6,i,k,j)
                c7=coeff(eqn,7,i,k,j)
                c4=coeff(eqn,4,i,k,j)
                c5=coeff(eqn,5,i,k,j)
                c8=coeff(eqn,8,i,k,j)
                c9=coeff(eqn,9,i,k,j)

                f1=vIN(i,k,j,eqn)
                f2=vIN(i-1,k,j,eqn)
                f3=vIN(i+1,k,j,eqn)
                f6=vIN(i,k,j-1,eqn)
                f7=vIN(i,k,jp1,eqn)
                f4=vIN(i,k-1,j,eqn)
                f5=vIN(i,k+1,j,eqn)
                gxy = (vIN(i+1,k,jp1,eqn)-vIN(i+1,k,j-1,eqn)&
                     -vIN(i-1,k,jp1,eqn)+vIN(i-1,k,j-1,eqn))
                gzy = (vIN(i,k+1,jp1,eqn)-vIN(i,k+1,j-1,eqn)&
                     -vIN(i,k-1,jp1,eqn)+vIN(i,k-1,j-1,eqn))

                addrate=zero
                do m=1,neqgas
                   addrate = addrate - Bm(eqn,m,i,k,j)*vIN(i,k,j,m)
                enddo
                rhsg = c2*f2+c3*f3+c6*f6+c7*f7 &
                     + c4*f4+c5*f5+c8*gxy+c9*gzy+addrate!&
                Vout(i,k,j,eqn) = -rhsg
             enddo
          enddo
       enddo
    enddo

    eqn =maxsize
    do j = drange(5)+1,drange(6)-1
       do k = drange(3),drange(4)
          do i = drange(1),drange(2)
             c2=csld(2,i,k,j)
             c3=csld(3,i,k,j)
             c6=csld(6,i,k,j)
             c7=csld(7,i,k,j)
             c4=csld(4,i,k,j)
             c5=csld(5,i,k,j)
             c8=csld(8,i,k,j)
             c9=csld(9,i,k,j)
             coe = csld(1,i,k,j)


             f2=vIN(i-1,k,j,eqn)
             f3=vIN(i+1,k,j,eqn)
             f6=vIN(i,k,j-1,eqn)
             f7=vIN(i,k,j+1,eqn)
             f4=vIN(i,k-1,j,eqn)
             f5=vIN(i,k+1,j,eqn)
             gxy = (+lambdas(i+1,k,j) * (vIN(i+1,k,j+1,eqn)&
                  - vIN(i+1,k,j-1,eqn))&
                  -lambdas(i-1,k,j) * (vIN(i-1,k,j+1,eqn)&
                  - vIN(i-1,k,j-1,eqn))&
                  +lambdas(i,k,j+1) * (vIN(i+1,k,j+1,eqn)&
                  - vIN(i-1,k,j+1,eqn))&
                  -lambdas(i,k,j-1) * (vIN(i+1,k,j-1,eqn)&
                  - vIN(i-1,k,j-1,eqn)))
             gzy = (+lambdas(i,k+1,j) * (vIN(i,k+1,j+1,eqn)&
                  - vIN(i,k+1,j-1,eqn))&
                  -lambdas(i,k-1,j) * (vIN(i,k-1,j+1,eqn)&
                  - vIN(i,k-1,j-1,eqn))&
                  +lambdas(i,k,j+1) * (vIN(i,k+1,j+1,eqn)&
                  - vIN(i,k-1,j+1,eqn))&
                  -lambdas(i,k,j-1) * (vIN(i,k+1,j-1,eqn)&
                  - vIN(i,k-1,j-1,eqn)))
             rhsg = c2*f2+c3*f3+c4*f4+c5*f5+c6*f6+c7*f7+c8*gxy+c9*gzy &
                  - coe*vIN(i,k,j,eqn) + csld(10,i,k,j)*vIN(i,k,0,maxsize)
             Vout(i,k,j,eqn) = -rhsg
          enddo
       enddo
    enddo

!boundary residuals
    Vout(:,:,drange(6),:) = 0d0
    Vout(:,:,drange(5),:) = 0d0
    if(implicitBC) then
       call BC_IMP(0,Vout,Vin,.false.)
    endif

    return
  end subroutine matvec
!----------------------------
  subroutine vecvec(vA,VB,Vprod)

    use global_data
    USE MYMPI

    implicit none

    real*8,Intent(IN) ::  vA(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize),VB(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize)
    real*8,Intent(OUT) :: Vprod
!local variables
    Integer :: i,k,M,j,it
    real*8 :: myvprod
!---------------------------------------

    myvprod = zero
    do M = 1,maxsize
       do j = drange(5),drange(6)
          do k = drange(3),drange(4)
             do i = drange(1),drange(2)
                myvprod = myvprod + VA(i,k,j,M)*VB(i,k,j,m)
             enddo
          enddo
       enddo
    enddo

    call mpi_allreduce(myVprod,Vprod,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)

    return
  end subroutine vecvec
!     **********************************************************************
!
!IMPLICIT BOUNDARY CONDITION SOLVER
!
!     **********************************************************************

  SUBROUTINE BC_IMP(color,Vres,Vfun,flag_I)

    USE GLOBAL_DATA
    USE BC_DATA
    USE MYMPI

    IMPLICIT NONE

!---------------------------------------------------------------------------
!  Dummy variables
    INTEGER, INTENT(IN) :: color
    real*8,  Intent(INOUT) :: Vres(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize)
    real*8,  Intent(INOUT) :: vfun(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize)
    logical,INTENT(IN),OPTIONAL :: flag_I
!  Local variables
    INTEGER :: i,k,j,iter,itmax,xs,ys,zs,xe,ye,ze,indx,eqn,eqA,eqZ
    REAL*8  :: dlg, const1, const2
    REAL*8  :: rmax,allrmax,prec_gss,ffold
    REAL*8  :: fx, fz
    REAL*8  :: biga,bigb,bigc,bigr,term
    LOGICAL :: flag
!---------------------------------------------------------------------------

    if(present(flag_I)) then
       flag = flag_I
    else
       flag = .true.
    end if

    xs = drange(1)
    xe = drange(2)
    ys = drange(5)
    ye = drange(6)
    zs = drange(3)
    ze = drange(4)

!...........................................................................
!  ZERO GRADIENT AT NORTH
!...........................................................................
    IF (color ==2 .and.flag) THEN
       do k = zs-1, ze+1
          do i = xs-1, xe+1
             do eqn = 1,maxsize-1
                Vfun(i,k,ny,eqn) = Vfun(i,k,ny-1,eqn)
             end do
          end do
       end do
    ENDIF
!
    eqn = 1
    j = 0
    if(flag) then
       TEMPERATURE: DO indx = 1 , finest_mesh%nclr(color)
          i =  finest_mesh%iclr(indx,color)
          k =  finest_mesh%kclr(indx,color)
          term = Vres(i,k,0,1) - capd(i,k)*Vfun(i,k,1,1) - capdd(i,k)*Vfun(i,k,2,1)  &
               -  cape(i,k)*Vfun(i,k,1,maxsize) - capee(i,k)*Vfun(i,k,2,maxsize)
          fx = Vfun(i+1,k,0,eqn)-Vfun(i-1,k,0,eqn)
          fz = Vfun(i,k+1,0,eqn)-Vfun(i,k-1,0,eqn)
          Vfun(i,k,0,eqn) = (capa(1,i,k)*fx + capb(1,i,k)*fz + term)/capc(1,i,k)
          Vfun(i,k,0,maxsize) = Vfun(i,k,0,eqn)
       ENDDO TEMPERATURE
    else
       do i = drange(1),drange(2)
          do k = drange(3),drange(4)
             term = - capd(i,k)*vfun(i,k,1,1) - capdd(i,k)*vfun(i,k,2,1)  &
                  -  cape(i,k)*vfun(i,k,1,maxsize) - capee(i,k)*vfun(i,k,2,maxsize)
             fx = vfun(i+1,k,0,eqn)-vfun(i-1,k,0,eqn)
             fz = vfun(i,k+1,0,eqn)-vfun(i,k-1,0,eqn)
             term =  term + capa(1,i,k)*fx + capb(1,i,k)*fz - capc(1,i,k)*vfun(i,k,0,eqn)
             Vres(i,k,j,1) = -term
             Vres(i,k,j,maxsize) = vfun(i,k,j,1) - vfun(i,k,j,maxsize)
          end do
       end do
    end if

!
!  SPECIES BOUNDARY CONDITIONS
! 
    eqA = 2
    eqZ = maxsize-1
    flagif: if(flag) then
       SPfun: DO EQN=eqA,eqZ
          DO indx = 1 , finest_mesh%nclr(color)
             i =  finest_mesh%iclr(indx,color)
             k =  finest_mesh%kclr(indx,color)
             term = Vres(i,k,0,eqn) - capd(i,k)*Vfun(i,k,1,eqn) - capdd(i,k)*Vfun(i,k,2,eqn) &
                  &                  + capT(eqn,i,k)*Vfun(i,k,0,1)
             fx=Vfun(i+1,k,0,eqn)-Vfun(i-1,k,0,eqn)
             fz=Vfun(i,k+1,0,eqn)-Vfun(i,k-1,0,eqn)
             Vfun(i,k,0,eqn) = (capa(2,i,k)*fx + capb(2,i,k)*fz + term)/capc(2,i,k)
          ENDDO
       ENDDO SPFUN
    else
       SPres: do EQN=eqA,eqZ
          do i = drange(1),drange(2)
             do k = drange(3),drange(4)
                term = - capd(i,k)*Vfun(i,k,1,eqn) - capdd(i,k)*Vfun(i,k,2,eqn) &
                     &                  + capT(eqn,i,k)*Vfun(i,k,0,1)
                fx=Vfun(i+1,k,0,eqn)-Vfun(i-1,k,0,eqn)
                fz=Vfun(i,k+1,0,eqn)-Vfun(i,k-1,0,eqn)
                term = term + capa(2,i,k)*fx + capb(2,i,k)*fz + term - capc(2,i,k)*vfun(i,k,0,eqn)
                Vres(i,k,j,EQN) = -term
             end do
          end do
       end do SPres
    endif flagif
!----------------------------------------------------------------------
    RETURN
  END SUBROUTINE BC_IMP
!**********************************************************************

end MODULE IMPLICIT
