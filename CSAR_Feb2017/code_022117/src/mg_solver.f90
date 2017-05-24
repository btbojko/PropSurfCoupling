MODULE mg_solver

  USE data_types
  USE global_data
  use timing

  IMPLICIT NONE

  PUBLIC  :: fmg_solve           ! solves the problem using the full multigrid method

  PRIVATE :: mg_cycle,        &  ! performs one complete multigrid cycle
       t_x_calc,        &  ! interpolates the correction from the coarse mesh to the fine mesh
       t_trans_r_calc      ! restricts the residual from the fine mesh to the coarse mesh

  character(LEN=1) :: BCTYPE
  INTEGER :: jmax

CONTAINS

  !***************************************************************************************************
  ! Subroutine to perform complete full multigrid cycles on the matrix equation Kx = f
  ! until convergence is achieved.  The process is described in the flowchart on page 45
  ! of: "Multigrid Methods: Fundamental Algorithms, Model Problem Analysis and 
  ! Applications" by K.STUBEN and U.TROTTENBERG, in "Lecture Notes in Mathematics 960:
  ! Multigrid  Methods", edited by W.Hackbusch and U.Trottenberg, Springer-Verlag 1982.
  !
  ! Note: the finest mesh is mesh #1, the coarsest mesh is mesh #NOMSHS.
  !***************************************************************************************************


!***************************************************************************************************

  SUBROUTINE fmg_solve(BCTYPE_IN)


!---------------------------------------------------------------------------------
!   Local variables:

    TYPE(mesh_pointers), POINTER :: current_mesh
    REAL(KIND=double) :: norm_f, norm_r, time

    character(LEN=1),OPTIONAL :: BCTYPE_IN

    INTEGER :: nxf,nzf,nyf,iter
!---------------------------------------------------------------------------------

    if(present(BCTYPE_IN)) then
       if(BCTYPE_IN == 'N' .or. BCTYPE_IN == 'D') then
          BCTYPE = BCTYPE_IN
       else
          BCTYPE = 'D'
       end if
    else
       BCTYPE = 'D'
    end if

    icycle = 0
!
!   Compute the load vector on each coarse mesh by restricting the load
!   

    current_mesh => finest_mesh
    converged = .FALSE.

    IF (iguess.eq.1 .OR. num_mg_meshes == 1) THEN

       do iter = 1, mg_maxcycle
          CALL gmres_solvep( current_mesh, nu_fine, .TRUE.)   !> USE .F. <
!!>          CALL gauss_relax( current_mesh, nu_fine, .TRUE.)   !> USE .F. <
          if(CONVERGED) EXIT
       enddo
       RETURN

    ELSEIF (iguess == 2) THEN

       CALL gmres_solvep( current_mesh, nu_fine, .FALSE.)   !> USE .F. <
!!>       CALL gauss_relax( current_mesh, nu_fine, .FALSE.)   !> USE .F. <

    ELSE


       DO

          IF ( .NOT. ASSOCIATED(current_mesh%coarse_mesh) ) EXIT

          nxf=current_mesh%xe +ngh;nzf=current_mesh%ze +ngh;nyf=current_mesh%ye +ngh
          wrk(-ngh:nxf,-ngh:nzf,-ngh:nyf) = current_mesh%f(-ngh:nxf,-ngh:nzf,-ngh:nyf)

!      CALL restrict( current_mesh%coarse_mesh )
          CALL transfer( current_mesh%coarse_mesh )

          current_mesh => current_mesh%coarse_mesh

       END DO

       CALL gmres_solvep( coarsest_mesh, nu1+nu2, .FALSE. )  !USE > .F. <
!!>       CALL gauss_relax( coarsest_mesh, nu1+nu2, .FALSE. )  !USE > .F. <

       !!    CALL output_results( coarsest_mesh )

       current_mesh => coarsest_mesh


       DO

          current_mesh => current_mesh%fine_mesh

          IF ( .NOT. ASSOCIATED(current_mesh) ) EXIT

!.................................................................................................
!     X_fine = zero
!     only up to nyf-1 assuming Dirichlet imposed at nyf and Neuman at 0
!     current_mesh%x(0:nxf,0:nzf,0:nyf-1) = zero      
!.................................................................................................
!
          nxf=current_mesh%xe +1;nzf=current_mesh%ze +1;nyf=current_mesh%ye +1
!      current_mesh%x(0:nxf,0:nzf,0:nyf-1) = zero      
          current_mesh%x(0:nxf,0:nzf,0:nyf) = zero      

!
!     Interpolate the solution on the coarse mesh to the fine mesh to give a good initial guess.
!
          CALL prolong( current_mesh%coarse_mesh, current_mesh )

          converged = .FALSE.

          CALL mg_cycle ( current_mesh )


       END DO    ! GUESS ON THE COARSEST MESH


    ENDIF   !iguess


    DO icycle = 1,mg_maxcycle-1

       IF ( converged ) EXIT

       current_mesh => finest_mesh  !START MESH

       CALL mg_cycle( current_mesh )

       IF ( converged ) EXIT

    END DO
    if(.not. converged) then
       CALL gmres_solvep( finest_mesh, nu_fine, .TRUE.) 
!!>    CALL gauss_relax( finest_mesh, nu_fine, .TRUE.) 
    endif
!--------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE fmg_solve
!********************************************************************************

  SUBROUTINE mg_cycle( start_mesh )

! Subroutine to perform one complete multigrid cycle using Jacobi
! preconditioned conjugate gradient relaxation.

!---------------------------------------------------------------------------------------------------

! Dummy arguments:

    TYPE(mesh_pointers), POINTER :: start_mesh

! Local variables:

    TYPE(mesh_pointers), POINTER :: current_mesh

    REAL(KIND=double) :: time

    INTEGER, POINTER, DIMENSION(:) :: c

    INTEGER :: start_mesh_number, current_mesh_number, error, i,icount

!---------------------------------------------------------------------------------------------------

    ALLOCATE( c(num_meshes), STAT=error ) 

! Perform one multigrid iteration step starting on mesh START_MESH.

    current_mesh => start_mesh

    start_mesh_number = start_mesh%mesh_num

    c = 0

    outer: DO

       DO 

          IF ( .NOT. ASSOCIATED(current_mesh%coarse_mesh) ) EXIT
          IF ( current_mesh%mesh_num == num_mg_meshes ) EXIT

          CALL gmres_solvep( current_mesh, nu1, .TRUE. )   !>USE .T.<
!!>          CALL gauss_relax( current_mesh, nu1, .TRUE. )   !>USE .T.<

          c(current_mesh%mesh_num) = c(current_mesh%mesh_num) + 1

          CALL restrict( current_mesh%coarse_mesh )

          current_mesh => current_mesh%coarse_mesh

          current_mesh%x = zero

       END DO

       CALL gmres_solvep( current_mesh, nu1+nu2, .FALSE. )  !>USE .F.<
!!>       CALL gauss_relax( current_mesh, nu1+nu2, .FALSE. )  !>USE .F.<

       current_mesh => current_mesh%fine_mesh

       IF ( .NOT. ASSOCIATED(current_mesh) ) RETURN

       DO

          CALL prolong( current_mesh%coarse_mesh, current_mesh )

          CALL gmres_solvep( current_mesh, nu2, .FALSE. )    !>USE .F.<
!!>          CALL gauss_relax( current_mesh, nu2, .FALSE. )    !>USE .F.<

          IF ( current_mesh%mesh_num == start_mesh_number ) EXIT outer

          IF ( c(current_mesh%mesh_num) /= gamma ) THEN
             CYCLE outer

          ELSE IF ( c(current_mesh%mesh_num) == gamma ) THEN
             c(current_mesh%mesh_num) = 0
          END IF

          current_mesh => current_mesh%fine_mesh

       END DO

    END DO outer

!---------------------------------------------------------------------------------------------------
  END SUBROUTINE mg_cycle
!***************************************************************************************************

  SUBROUTINE t_x_calc( coarse_mesh, fine_mesh ) 

! Subroutine to interpolate the correction from the coarse mesh to
! the fine mesh by computing the product [T]{x} using an element-by-
! element approach on the coarse mesh.

! The coarse mesh correction is stored in COARSE_MESH%X
! The new fine mesh solution is stored in FINE_MESH%X

!---------------------------------------------------------------------------------------------------

! Dummy arguments:

    TYPE(mesh_pointers), POINTER :: coarse_mesh
    TYPE(mesh_pointers), POINTER :: fine_mesh

! Local variables:

    REAL(KIND=double), DIMENSION(:), POINTER :: x_fine_temp

    INTEGER :: error

!---------------------------------------------------------------------------------------------------

! Assign and zero the array that will contain the interpolated fine mesh correction.

    !!x_fine_temp => work_vectors%kp

!  x_fine_temp(1:fine_mesh%num_equ_free) = zero

! Compute the contribution to [T]{x} from each of the element groups.

    write(*,*)'INTERPOLATING FROM COARSE MESH TO FINE, GOING UP',&
         coarse_mesh%mesh_num,fine_mesh%mesh_num,coarse_mesh%xe,fine_mesh%xe,&
         coarse_mesh%ye,fine_mesh%ye
    !!IF ( ASSOCIATED(coarse_mesh%b8_bbar) ) CALL b8_bbar_t_x_calc( coarse_mesh )
    !!IF ( ASSOCIATED(coarse_mesh%b8_me)   ) CALL b8_me_t_x_calc( coarse_mesh )

! Add the interpolated coarse mesh correction to the current fine mesh solution.

    !!fine_mesh%x(1:fine_mesh%num_equ_free) = fine_mesh%x(1:fine_mesh%num_equ_free) +   &
    !!                                        x_fine_temp(1:fine_mesh%num_equ_free)

    NULLIFY (x_fine_temp)

!---------------------------------------------------------------------------------------------------

  END SUBROUTINE t_x_calc

!***************************************************************************************************

  SUBROUTINE t_trans_r_calc( coarse_mesh, time )

! Subroutine to restrict the residual from the fine mesh to the
! coarse mesh by computing the product [T]^T{r} using an element-by-
! element approach on the coarse mesh.

! The fine   mesh residual is stored in WORK_VECTORS%R
! The coarse mesh residual is stored in COARSE_MESH%F

!---------------------------------------------------------------------------------------------------

! Dummy arguments:

    TYPE(mesh_pointers), POINTER :: coarse_mesh

    REAL(KIND=double) :: time

!---------------------------------------------------------------------------------------------------

! Initialize the residual on the coarse mesh.

!  coarse_mesh%f(1:coarse_mesh%num_equ_free) = zero

! Compute the contribution to [T]^T{r} from each of the element groups.

    write(*,*)'RESTRICTING FROM FINE MESH TO COARSE, GOING DOWN',&
         coarse_mesh%mesh_num,&
         coarse_mesh%xe,coarse_mesh%ye
!    IF ( ASSOCIATED(coarse_mesh%b8_bbar) ) CALL b8_bbar_t_trans_r_calc( coarse_mesh, time )
!    IF ( ASSOCIATED(coarse_mesh%b8_me)   ) CALL b8_me_t_trans_r_calc( coarse_mesh, time )

!---------------------------------------------------------------------------------------------------

  END SUBROUTINE t_trans_r_calc
!*********************************************************************************

  SUBROUTINE gauss_relax ( current_mesh, nu, resid_flag )

    USE GLOBAL_DATA 
    USE MYMPI

    IMPLICIT NONE

    include "parallel.h"

! Subroutine performs NU cycles of GS
! relaxation on mesh CURRENT_MESH.
! NOTE ASSUME NO PROCESSOR_PARTITIONING IN Y
! a 2 colors scheme works  if there is no partitioning in y direction
! there is no cross term in the x-z plane. This implementation also requires
! that the  number of pts in x and z is even
! This is necessarily true except on the coarsest grid

!---------------------------------------------------------------------------------------------------
!   Dummy variables:

    TYPE(mesh_pointers) :: current_mesh
!    TYPE(mesh_pointers), POINTER :: current_mesh

    INTEGER, INTENT(IN) :: nu                ! number of relaxations 
    LOGICAL, INTENT(IN) :: resid_flag        ! if this is true it evaluates the residual

!   Local Variables
    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: px
    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: fx

    type(cppevec), DIMENSION(:,:), POINTER :: CPPE
    REAL(KIND=double), DIMENSION(:,:), POINTER :: ccmg
    REAL(KIND=double), DIMENSION(:), POINTER :: ey

    INTEGER :: icount,i,j,k,xe,ze,ye,nuu
    INTEGER :: indx,icolor,rcolor, tmpcomm
    INTEGER :: js, jrflo_converto_p,jp1

    REAL(KIND=double) :: resmax,myres,omg_ppe,cc(12),vv(12)
    CHARACTER(LEN=8) :: msg
!---------------------------------------------------------------------------------------------------

!
!   RETURN immediatly if the processor has no points
!
    if(current_mesh%blank) RETURN

    px => current_mesh%x ; fx => current_mesh%ff
    CPPE    =>  current_mesh%cppe

    xe = current_mesh%xe ; ze = current_mesh%ze; ye = current_mesh%ye 
    tmpcomm = current_mesh%comm3d
    nuu=nu; if(current_mesh%mesh_num <= 1 .and. nu > 0) nuu=nu_fine
    omg_ppe=1.15d0

!
!   START ITERATIONS
!
    jmax = current_mesh%ye+1

    countdo: do  icount = 1,nuu
       

       colordo: do icolor=1,2


          rcolor=mod(icolor+2,2)+1     !2,1

          call color_swap(current_mesh,rcolor)
!!>

          indxdo: do indx = 1 , current_mesh%nclr(icolor)


             i = current_mesh%iclr(indx,icolor)
             k = current_mesh%kclr(indx,icolor)
             ccmg => cPPE(i,k)%ccmg

             do j=js,jmax-1
                if(BCTYPE == 'D'.or. j< jmax-1) then
                   jp1 = j+1
                else
                   jp1 = j
                end if
!!>                jp1 = merge(j+1,j,BCTYPE == 'D'.or. j< jmax-1)


                fy(j) = fx(i,k,j) + ccmg(2,j)*px(i-1,k,j)+                     &
                     ccmg(3,j)*px(i+1,k,j)  +  ccmg(4,j)*px(i,k-1,j)+          &
                     ccmg(5,j)*px(i,k+1,j)  +  ccmg(8,j)*px(i-1,k,j-1) +       &
                     ccmg(9,j)*px(i+1,k,j-1)+  ccmg(10,j)*px(i+1,k,jp1)+       &
                     ccmg(11,j)*px(i-1,k,jp1)+ ccmg(12,j)*px(i,k-1,j-1)+       &
                     ccmg(13,j)*px(i,k+1,j-1)+ ccmg(14,j)*px(i,k+1,jp1)+       &
                     ccmg(15,j)*px(i,k-1,jp1)

                ay(j) = ccmg(6,j)
                cy(j) = ccmg(7,j)
                by(j) = ccmg(1,j)

             enddo
             if(BCTYPE == 'D') then
                fy(jmax) = 0.0d0
                ay(jmax) = 0.0d0
                cy(jmax) = 0.0d0
                by(jmax) = 1.0d0
             elseif(BCTYPE == 'N') then
                fy(jmax) = 0.0d0
                ay(jmax) = -1.0d0
                cy(jmax) = 0.0d0
                by(jmax) =  1.0d0
             end if

! --- THOMAS algorithm
             do j=js+1,jmax
                ay(j) = ay(j)/by(j-1)
                by(j) = by(j)-ay(j)*cy(j-1)
                fy(j) = fy(j)-ay(j)*fy(j-1)
             enddo

             fy(jmax) = fy(jmax)/by(jmax)
             px(i,k,jmax) = fy(jmax)

             do j=jmax-1,js,-1
                fy(j) = (fy(j)-cy(j)*fy(j+1))/by(j)
                px(i,k,j) = fy(j)
             enddo
             
          enddo indxdo   !do indx
       enddo colordo   ! do icolor
    enddo countdo   !END   do  icount = 1,nuu

    if(BCTYPE == 'N') then
       px(:,:,ye+1) = px(:,:,ye)
    end if

!
    rcolor=2
    call color_swap(current_mesh,rcolor)

    IF(resid_flag) then
!
!EVAL residual for transfering to lower levels
!
       myres = zero

       do i=0,xe
          do k=0,ze

             ccmg => cPPE(i,k)%ccmg
             do j=0,jmax-1
                jp1 = merge(j+1,j,BCTYPE == 'D'.or. j< jmax-1)


                wrk(i,k,j) = fx(i,k,j) -  ccmg(1,j)*PX(i,k,j) + &
                  & ccmg(2,j)*PX(i-1,k,j) + ccmg(3,j)*PX(i+1,k,j) + ccmg(4,j)*PX(i,k-1,j)&
                  &+ ccmg(5,j)*PX(i,k+1,j) - ccmg(6,j)*PX(i,k,j-1) - ccmg(7,j)*PX(i,k,jp1) + ccmg(8,j)*PX(i-1,k,j-1)&
                  &+ ccmg(9,j)*PX(i+1,k,j-1)+ ccmg(10,j)*PX(i+1,k,jp1) + ccmg(11,j)*PX(i-1,k,jp1) + ccmg(12,j)*PX(i,k-1,j-1)&
                  &+ ccmg(13,j)*PX(i,k+1,j-1)+ ccmg(14,j)*PX(i,k+1,jp1)+ ccmg(15,j)*PX(i,k-1,jp1)

                myres = max (myres, abs(wrk(i,k,j)))

                if(wrk(i,k,j) > 1d40) then
!!>                if(myid == 0.and.current_mesh%mesh_num ==6) then
!!>                   print'(4i3,1p2e12.4,1(a,12e12.4),3e12.4,1(a,12e12.4))',icycle,i,k,j,wrk(i,k,j),fx(i,k,j),'CC',cc, cPPE(i,k)%ccmg(6:7,j),cPPE(i,k)%ccmg(1,j),'vv',vv
                endif

             enddo
             if(BCTYPE == 'D') then
                wrk(i,k,jmax) = zero
             elseif(BCTYPE == 'N') then
                wrk(i,k,jmax) = wrk(i,k,jmax-1)
             end if
          enddo
!!>          if(myid == 0)print '(i3,1p9999e11.3)',i,wrk(i,0,0:jmax-1)
       enddo
       if(current_mesh%mesh_num == 1) then
          msg(1:8) ='<<<<<<<<<<<<<<<=1'
       else
          msg(1:8) ='###############>1'
       endif
       if(current_mesh%mesh_num == 1) then
          CALL MPI_ALLREDUCE(myres,resmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,tmpcomm,ierr)
          converged = resmax < prec_ppe
          if(myid==0 .and. ipack > -1 .and. mod(ncyc,writemod) == 0) print '(i3,1pe12.4,1x,L1,i3,a)',&
               &icycle,resmax,converged,current_mesh%mesh_num,TRIM(msg)
          maxres(3) = resmax
       endif

    ENDIF   !if (resid_flag)

!---------------------------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE gauss_relax
!***************************************************************************************************


!***************************************************************************************************
  SUBROUTINE restrict ( coarse_mesh )

!fine_mesh -->> coarse mesh     T_trans_r_calc (ALI)

    USE GLOBAL_DATA 
    USE MYMPI

    IMPLICIT NONE

    include "parallel.h"

!---------------------------------------------------------------------------------------------------
! Dummy variables:

    TYPE(mesh_pointers) :: coarse_mesh

! Local Variables

    TYPE(mesh_pointers), POINTER :: fine_mesh

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: px

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: fx

    INTEGER :: i,j,k,id,kd
    INTEGER :: ic,kc,jc,mn, col(6)
    INTEGER :: ncx,ncy,ncz,xe,ze,ye

    REAL(KIND=double) :: rm,rj,rp
!---------------------------------------------------------------------------------------------------
!
! This implementation uses full weghted residual restriction which is (little) 
! different from what is used in rocsolid
!
    fine_mesh => coarse_mesh%fine_mesh
    if (fine_mesh%blank) RETURN

    px =>  coarse_mesh%x ; fx =>  coarse_mesh%f
    ncx = coarse_mesh%xe + 1; ncy = coarse_mesh%ye + 1 ;ncz = coarse_mesh%ze + 1;
    xe = fine_mesh%xe ; ye = fine_mesh%ye; ze = fine_mesh%ze;
    mn = fine_mesh%mesh_num; id = id_mg(mn); kd = kd_mg(mn)


    col = 0;col(2) = xe;col(4) = ze; col(6) = ye;

    if ( iscoarse_xz(mn) ) &
         call parallel_swap(wrk,col, fine_mesh)

    if(coarse_mesh%blank) RETURN


    if ( iscoarse_xz(mn) .AND. iscoarse_y(mn) ) then
!
!   coarsening in x,y,z
!
       do kc=0,ncz-1
          k = kc+kc+kd
          do ic=0,ncx-1
             i = ic+ic+id

!
!       weight on j-1,j,j+1 z planes in rm,rk,rp
!
             DO j=0,ye
                slambda(j) = ( wrk(i-1,k-1,j)+ wrk(i+1,k-1,j)+ wrk(i-1,k+1,j)+&
                     wrk(i+1,k+1,j)+2.*( wrk(i-1,k,j)+ wrk(i+1,k,j)+&
                     wrk(i,k-1,j)+ wrk(i,k+1,j))+4.* wrk(i,k,j))*.0625
             ENDDO

             fx(ic,kc,0) = third*(slambda(1) + two*slambda(0) )

             DO jc = 1,ncy-1
                j = jc+jc
                fx(ic,kc,jc) = 0.25*(slambda(j-1) + 2.0d0* slambda(j) + slambda(j+1))
             ENDDO

          end do
       end do


    elseif ( iscoarse_xz(mn) ) then

!
!   coarsening in x,z but not y, 2D coasening
!
       do kc=0,ncz-1
          k = kc+kc+kd
          do jc=0,ncy-1
             j = jc
             do ic=0,ncx-1
                i = ic+ic+id

                rj=( wrk(i-1,k-1,j)+ wrk(i+1,k-1,j)+ wrk(i-1,k+1,j)+&
                     wrk(i+1,k+1,j)+2.*( wrk(i-1,k,j)+ wrk(i+1,k,j)+&
                     wrk(i,k-1,j)+ wrk(i,k+1,j))+4.* wrk(i,k,j))*.0625

                fx(ic,kc,jc) = rj

             end do
          end do
       end do

    elseif ( iscoarse_y(mn) ) then
!
!   coarsening in y but not xy, 1D coasening
!
       do kc=0,ncz-1
          k = kc
          do ic=0,ncx-1
             i = ic
             do jc=1,ncy-1
                j = jc+jc

                fx(ic,kc,jc) = quarter*(wrk(i,k,j+1)+two*wrk(i,k,j)+wrk(i,k,j-1)) 

             end do
             fx(ic,kc,0) =  third*(wrk(i,k,1)+two*wrk(i,k,0))
          end do
       end do


    else

       write(*,*) 'NO COARSENING, CHECK THE CALL!!'

    endif

!----------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE restrict
!**********************************************************************************


!**********************************************************************************
  SUBROUTINE prolong (coarse_mesh ,fine_mesh)

    USE GLOBAL_DATA 
    USE MYMPI

    IMPLICIT NONE

    include "parallel.h"
!---------------------------------------------------------------------------------------------------
! Dummy variables:

    TYPE(mesh_pointers) :: coarse_mesh,fine_mesh


! Local Variables

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: cx, fx    !fx is the solution on the fine mesh

    INTEGER :: i,j,k
    INTEGER :: ic,kc,jc,mn,id,kd
    INTEGER :: ncx,ncy,ncz
    INTEGER :: nfx,nfy,nfz

    REAL(KIND=double) :: rm,rj,rp
!---------------------------------------------------------------------------------------------------

    if(fine_mesh%blank) RETURN
!
!   initialize pointers and LOOP limits
!
    cx => coarse_mesh%x ; fx => fine_mesh%x  
    ncx = coarse_mesh%xe + 1; ncy = coarse_mesh%ye + 1 ;ncz = coarse_mesh%ze + 1;
    nfx = fine_mesh%xe + 1; nfy = fine_mesh%ye + 1 ;nfz = fine_mesh%ze + 1;
    mn = fine_mesh%mesh_num; id = id_mg(mn); kd = kd_mg(mn)

!
! ASSUMING pdims(2) > pdims(1) and nx=nz these two are equivalent
!
!    if(any(coarse_mesh%pskip) .AND. iscoarse_xz(mn)) &
    if(coarse_mesh%pskip(2) .AND. iscoarse_xz(mn)) &
         CALL fill_voids(coarse_mesh,fine_mesh)

    if ( iscoarse_xz(mn) .AND. iscoarse_y(mn) ) then
!
!   Interpolating in x,y,z
!   It is done in three step, 1-D interp, 2-D interp and finally 3-D interp
!
       do kc=-kd,ncz
          k = kc+kc+kd
          do jc=0,ncy
             j = jc+jc
             do ic=-id,ncx
                i = ic+ic+id

                wrk(i,k,j) = cx(ic,kc,jc)
                if(ic.lt.ncx)&
                     wrk(i+1,k,j) = 0.5*(cx(ic,kc,jc)+cx(ic+1,kc,jc))
                if(kc.lt.ncz)&
                     wrk(i,k+1,j) = 0.5*(cx(ic,kc,jc)+cx(ic,kc+1,jc))
                if(jc.lt.ncy)&
                     wrk(i,k,j+1) = 0.5*(cx(ic,kc,jc)+cx(ic,kc,jc+1))

             enddo
          enddo
       enddo


       do kc=-kd,ncz
          k = kc+kc+kd
          do jc=0,ncy
             j = jc+jc
             do ic=-id,ncx
                i = ic+ic+id

                if(ic.lt.ncx.and.kc.lt.ncz)&
                     wrk(i+1,k+1,j) = 0.5*(wrk(i+1,k,j)+wrk(i+1,k+2,j))
                if(ic.lt.ncx.and.jc.lt.ncy)&
                     wrk(i+1,k,j+1) = 0.5*(wrk(i+1,k,j)+wrk(i+1,k,j+2))
                if(kc.lt.ncz.and.jc.lt.ncy)&
                     wrk(i,k+1,j+1) = 0.5*(wrk(i,k+1,j)+wrk(i,k+1,j+2))

             enddo
          enddo
       enddo

       do kc=-kd,ncz-1
          k = kc+kc+kd
          do jc=0,ncy-1
             j = jc+jc
             do ic=-id,ncx-1
                i = ic+ic+id

                wrk(i+1,k+1,j+1) = 0.5*(wrk(i+1,k,j+1)+wrk(i+1,k+2,j+1))

             enddo
          enddo
       enddo

    elseif( iscoarse_xz(mn)) then
!
!   interpolating in x,z
!   It is done in two steps, 1-D interp, 2-D interp 
!
       do kc=-kd,ncz
          k = kc+kc+kd
          do jc=0,ncy
             j = jc
             do ic=-id,ncx
                i = ic+ic+id

                wrk(i,k,j) = cx(ic,kc,jc)
                if(ic.lt.ncx)&
                     wrk(i+1,k,j) = 0.5*(cx(ic,kc,jc)+cx(ic+1,kc,jc))
                if(kc.lt.ncz)&
                     wrk(i,k+1,j) = 0.5*(cx(ic,kc,jc)+cx(ic,kc+1,jc))
! 
             enddo
          enddo
       enddo


       do jc=0,ncy
          j = jc
          do kc=-kd,ncz-1
             k = kc+kc+kd
             do ic=-id,ncx-1
                i = ic+ic+id

                wrk(i+1,k+1,j) = 0.5*(wrk(i+1,k,j)+wrk(i+1,k+2,j))

             enddo
          enddo
       enddo
!wrk(2*ncx-1,2*ncz-1,:)= 0.5*(cx(ncx,ncz-1,:)+cx(ncx-1,ncz,:))

    ELSEIF ( iscoarse_y(mn) ) then
!
!   interpolating in y only
!
       do kc=0,ncz
          k = kc
          do ic=0,ncx
             i = ic
             do jc=0,ncy-1
                j = jc+jc

                wrk(i,k,j) = cx(ic,kc,jc)
                wrk(i,k,j+1) = 0.5*(cx(ic,kc,jc)+cx(ic,kc,jc+1))

             enddo
          enddo
       enddo

    ELSE

       write(*,*) 'NO COARSENING, CHECK THE CALL!!'

    ENDIF

!
!   Inject the coarse grid correction into the fine mesh approximation
!
    do i=0,nfx-1
       do k=0,nfz-1
          do j=0,nfy-1

             fx(i,k,j) = fx(i,k,j)+wrk(i,k,j)

          enddo
       enddo
    enddo

!---------------------------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE prolong
!***************************************************************************************************


!***************************************************************************************************
  SUBROUTINE transfer ( coarse_mesh )

!This subroutine is the same as restrict (more or less). It is separate (for now)
! because at the start of the algorithm we just transfer the known term down

    USE GLOBAL_DATA 

    IMPLICIT NONE

!---------------------------------------------------------------------------------------------------
! Dummy variables:

    TYPE(mesh_pointers) :: coarse_mesh

! Local Variables

    TYPE(mesh_pointers), POINTER :: fine_mesh

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: cx
    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: cf
    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: fx

    INTEGER :: i,j,k,id,kd
    INTEGER :: ic,kc,jc,mn
    INTEGER :: ncx,ncy,ncz,incy,incxz
    INTEGER :: xe,ze,ye
!---------------------------------------------------------------------------------------------------
!
!   RETURN IMMEDIATELY IF THE PROC HAS NO POINTS
!
    if(coarse_mesh%blank) RETURN

    fine_mesh => coarse_mesh%fine_mesh
    cx => coarse_mesh%x ; cf =>  coarse_mesh%f
    fx =>   fine_mesh%x
    ncx = coarse_mesh%xe + 1; ncy = coarse_mesh%ye + 1 ;ncz = coarse_mesh%ze + 1;
    mn = fine_mesh%mesh_num; id = id_mg(mn); kd = kd_mg(mn)

!
    if ( iscoarse_xz(mn) .AND. iscoarse_y(mn) ) then
       incxz=2
       incy=2
    else if (iscoarse_xz(mn)) then
       incxz=2
       incy=1
    else if (iscoarse_y(mn)) then
       incxz=1
       id=0;kd=0
       incy=2
    else
       write(*,*) 'ERROR in TRANSFER, mesh # = ',mn
    endif
!
!   coarsening in x,y,z
!
    do kc=0,ncz-1
       k = kc*incxz+kd
       do ic=0,ncx-1
          i = ic*incxz+id
          do jc=0,ncy-1
             j = jc*incy
!
! JUST inject
!
             cf(ic,kc,jc) = wrk(i,k,j)
             cx(ic,kc,jc) =  fx(i,k,j)

          end do
       end do
    end do

!
!   This extra communication is necessary only on the corasest grid
!   when iex is odd
!
    if (mn == num_mg_meshes-1 .AND. isoddcoarse)&
         call color_swap(coarse_mesh,1)

!---------------------------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE transfer
!***************************************************************************************************


!***************************************************************************************************
  SUBROUTINE output_results( current_mesh )

    USE GLOBAL_DATA 

    IMPLICIT NONE

!---------------------------------------------------------------------------------------------------
! Dummy variables:

    TYPE(mesh_pointers) :: current_mesh


! Local Variables

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: px

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: fx

    INTEGER :: i,j,k
    INTEGER :: ic,kc,jc,mn
    INTEGER :: ncx,ncy,ncz,incy
!---------------------------------------------------------------------------------------------------

    px =>  current_mesh%x ; fx =>  current_mesh%f
    ncx = current_mesh%xe + 1; ncy = current_mesh%ye + 1 ;ncz = current_mesh%ze + 1;
    mn = current_mesh%mesh_num

    do jc=0,ncy-1
       do kc=0,ncz-1
          do ic=0,ncx-1

             if(myid == 0) &
                  write(*,*)myid,ic,kc,jc,'OUT',px(ic,kc,jc),fx(ic,kc,jc),ncy

          end do
       end do
    end do

!---------------------------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE output_results
!***************************************************************************************************

!!>
!***************************************************************************************************
  subroutine GMres_solveP(current_mesh, nu, resid_flag)
!
    use global_data
    USE MYMPI
    use Gather_R

    implicit none
    include 'parallel.h'

    TYPE(mesh_pointers) :: current_mesh
    Integer,Intent(IN) :: nu
    logical,Intent(IN) :: resid_flag
!local variables
    Integer :: i,k,j,M,it,ii,kk,jj,col(6),nugmres,mn,maxGMREScyc
    real*8 :: betaGMRES,alphaGMRES,htmp,rd,hd,rr,initialresid,TOLGMRES
    logical :: stoppage
    character(len =5) :: prec
!---------------------------------------

    col(5)=drange(5);col(6)=drange(6);
    col(3)=1;col(4)=maxsize
    col(1) = 1
    col(2) = 2
    prec = 'left'
    nugmres = 1
    TOLGMRES = 1d-7
    mn = current_mesh%mesh_num
    maxGMREScyc = merge(2000,15,mn==num_mg_meshes)

    if(current_mesh%blank) RETURN

    call color_swap(current_mesh,1)
    call color_swap(current_mesh,2)
    do j = 0,current_mesh%ye
       do k = -1,current_mesh%ze+1
          do i =-1,current_mesh%xe+1
             current_mesh%xold(i,k,j) = current_mesh%x(i,k,j)
             current_mesh%gm(1)%mv(i,k,j) = current_mesh%x(i,k,j)
          end do
       end do
    end do

!!>    if(myid == 0) print*,icycle,maxval(abs(current_mesh%x))
!!>    if(mn == 1) call gatherbuffMG(current_mesh,current_mesh%x,'XSOLMG',icycle)

    call matvecp(current_mesh,current_mesh%gm(1)%mv,current_mesh%ff)

    do j = 0,current_mesh%ye
       do k = -1,current_mesh%ze+1
          do i =-1,current_mesh%xe+1
             current_mesh%ff(i,k,j) = current_mesh%f(i,k,j) - current_mesh%ff(i,k,j)
          end do
       end do
    end do

!!>    if(mn == 1) call gatherbuffMG(current_mesh,current_mesh%ff,'XRESMG',icycle)
    call vecvecp(current_mesh,current_mesh%ff,current_mesh%ff,betaGMRES)
    betaGMRES = max(sqrt(betaGMRES),1d-12)
!!>    if(myid == 0) print*,'betaGMRESINITIAL',betaGMRES,mn
    
    if(prec == 'left ') then
       current_mesh%x = 0d0
       call gauss_relax(current_mesh,nugmres,.FALSE.)   !p initial guess should be =0
    else
       do j = 0,current_mesh%ye
          do k = -1,current_mesh%ze+1
             do i =-1,current_mesh%xe+1
                current_mesh%x(i,k,j) = current_mesh%ff(i,k,j)
             end do
          end do
       end do
       call color_swap(current_mesh,1)
       call color_swap(current_mesh,2)
    end if
       
!
    call vecvecp(current_mesh,current_mesh%x,current_mesh%x,betaGMRES)
    betaGMRES = max(sqrt(betaGMRES),1d-12)
    
    do j = 0,current_mesh%ye
       do k = -1,current_mesh%ze+1
          do i =-1,current_mesh%xe+1
             current_mesh%gm(1)%mv(i,k,j) = current_mesh%x(i,k,j)/betaGMRES
          enddo
       enddo
    enddo
    dgm(1) = betaGMRES;
    initialresid = betaGMRES
!
    call vecvecp(current_mesh,current_mesh%f,current_mesh%f,alphaGMRES)
    alphaGMRES = max(sqrt(alphaGMRES),1d-12)


    kk = 0
    GmresLOOP:do while (.true.)
    
       kk = kk+1
       stoppage = kk >= min(ubound(current_mesh%gm,1),maxGMREScyc)-1

       call matvecp(current_mesh,current_mesh%gm(kk)%mv,current_mesh%ff)   !it needs ghost points but not BC points

       if(prec == 'left ') then
          current_mesh%x = zero          
          call gauss_relax(current_mesh,nugmres,.FALSE.)
       else
          do j = 0,current_mesh%ye
             do k = -1,current_mesh%ze+1
                do i =-1,current_mesh%xe+1
                   current_mesh%x(i,k,j) = current_mesh%ff(i,k,j)
                end do
             end do
          end do
          call color_swap(current_mesh,1)
          call color_swap(current_mesh,2)
       end if

       do ii = 1,kk
          call vecvecp(current_mesh,current_mesh%x,current_mesh%gm(ii)%mv,hgm(ii,kk))
       enddo

       m = 1
       do ii = 1,kk
          do j = 0,current_mesh%ye
             do k = -1,current_mesh%ze+1
                do i = -1,current_mesh%xe+1
                   current_mesh%x(i,k,j) = current_mesh%x(i,k,j) - hgm(ii,kk)*current_mesh%gm(ii)%mv(i,k,j)
                enddo
             enddo
          enddo
       enddo

       call vecvecp(current_mesh,current_mesh%x,current_mesh%x,betaGMRES)
       betaGMRES = max(sqrt(betaGMRES),1D-12)
       hgm(kk+1,kk) = betaGMRES

       do j =0,current_mesh%ye
          do k = -1,current_mesh%ze+1
             do i =  -1,current_mesh%xe+1
                current_mesh%gm(kk+1)%mv(i,k,j) = current_mesh%x(i,k,j)/betaGMRES
             enddo
          enddo
       enddo
       
       do ii = 1,kk-1
          htmp = c1gm(ii) * hgm(ii,kk) + c2gm(ii) * hgm(ii+1,kk);
          hgm(ii+1,kk) = - c2gm(ii) * hgm(ii,kk) + c1gm(ii) * hgm(ii+1,kk);
          hgm(ii,kk) = htmp;
       end do
       rD = hgm(kk,kk);
       hD = hgm(kk+1,kk);
       rr = max(sqrt(rD*rD+hD*hD),1D-12)
       c1gm(kk) = rD/rr;
       c2gm(kk) = hD/rr;
       hgm(kk,kk) = rr;
       hgm(kk+1,kk) = zero;
       dgm(kk+1) = -c2gm(kk)*dgm(kk);
       dgm(kk) = c1gm(kk)*dgm(kk);
       errGMres = abs(dgm(kk+1))!/initialresid
!!>       if(myid == 0 .and. mn == 2)write(*,*)icycle,kk+1,'GMresAbs',abs(dgm(kk+1))
       if (errGMres  < TOLGMRES.or. stoppage) then
          do ii = kk,1,-1
             ygm(ii) = dgm(ii)/hgm(ii,ii);
             do jj = ii-1,1,-1
                dgm(jj) = dgm(jj) - hgm(jj,ii)*ygm(ii);
             end do
             do j = 0,current_mesh%ye
                do k = -1,current_mesh%ze+1
                   do i = -1,current_mesh%xe+1
                      current_mesh%xold(i,k,j) = current_mesh%xold(i,k,j) + ygm(ii)*current_mesh%gm(ii)%mv(i,k,j) 
                   enddo
                enddo
             enddo
          end do
          exit GmresLOOP !done
       endif
    enddo GmresLOOP

    do j = 0,current_mesh%ye
       do k = -1,current_mesh%ze+1
          do i =-1,current_mesh%xe+1
             current_mesh%x(i,k,j) = current_mesh%xold(i,k,j)
             current_mesh%ff(i,k,j) = current_mesh%f(i,k,j)
          end do
       end do
    end do
    call color_swap(current_mesh,1)
    call color_swap(current_mesh,2)
    call gauss_relax(current_mesh,nugmres,resid_flag)



!!>    do j = 0,current_mesh%ye
!!>       do k = -1,current_mesh%ze+1
!!>          do i =-1,current_mesh%xe+1
!!>             current_mesh%xold(i,k,j) = current_mesh%x(i,k,j)
!!>             current_mesh%mv(i,k,j,1) = current_mesh%x(i,k,j)
!!>          end do
!!>       end do
!!>    end do
!!>    call matvecp(current_mesh,current_mesh%mv,current_mesh%ff,1)
!!>    do j = 0,current_mesh%ye
!!>       do k = -1,current_mesh%ze+1
!!>          do i =-1,current_mesh%xe+1
!!>             current_mesh%ff(i,k,j) = current_mesh%f(i,k,j) - current_mesh%ff(i,k,j)
!!>          end do
!!>       end do
!!>    end do
!!>    call vecvecp(current_mesh,current_mesh%ff,current_mesh%ff,betaGMRES)
!!>    betaGMRES = max(sqrt(betaGMRES),1d-12)

!!>    if(myid == 0) print*,'betaGMRESFINAL',betaGMRES,mn

!!>    if(mn == 1.and.betaGMRES<=800000.d0) call gatherbuffMG(current_mesh,current_mesh%x,'XFINMG',icycle)
!!>    if(betaGMRES<=800000.d0) call gatherbuffMG(current_mesh,current_mesh%x,'XFINMG',icycle+mn*10)


    return
  end subroutine GMres_solveP
!----------------------------

!----------------------------
  subroutine matvecp(current_mesh,vIN,vOUT)

    use global_data
    USE MYMPI 
    USE BC_DATA

    implicit none

    TYPE(mesh_pointers) :: current_mesh
    real*8,POINTER ::  vIN (:,:,:)
    real*8,POINTER :: vOUT(:,:,:)
!local variables
    Integer :: i,k,j,M,it,jp1,eqn
    Real*8 :: tmp(maxsize),mynormresid
    Real*8 :: c2,c3,c6,c7,c4,c5,c8,c9,f2,f1,f3,f6,f7,f4,f5,gxy 
    Real*8 :: gzy,addrate,rr,coe,rhsg,fx,fz,term
    REAL(KIND=double) :: resmax,myres,omg_ppe,cc(12),vv(12)
    real(KIND=double),POINTER :: ccmg(:,:)
!---------------------------------------

    jmax = current_mesh%ye+1


    do i=0,current_mesh%xe
       do k=0,current_mesh%ze

          ccmg => current_mesh%cPPE(i,k)%ccmg
          do j=0,current_mesh%ye
             jp1 = merge(j+1,j,BCTYPE == 'D'.or. j< jmax-1)

             vOUT(i,k,j)  =  ccmg(1,j)*VIN(i,k,j) - ccmg(2,j)*VIN(i-1,k,j) - ccmg(3,j)*VIN(i+1,k,j) - ccmg(4,j)*VIN(i,k-1,j)&
                  &- ccmg(5,j)*VIN(i,k+1,j) + ccmg(6,j)*VIN(i,k,j-1) + ccmg(7,j)*VIN(i,k,jp1) - ccmg(8,j)*VIN(i-1,k,j-1)&
                  &- ccmg(9,j)*VIN(i+1,k,j-1)- ccmg(10,j)*VIN(i+1,k,jp1) - ccmg(11,j)*VIN(i-1,k,jp1) - ccmg(12,j)*VIN(i,k-1,j-1)&
                  &- ccmg(13,j)*VIN(i,k+1,j-1)- ccmg(14,j)*VIN(i,k+1,jp1)- ccmg(15,j)*VIN(i,k-1,jp1)

          enddo
       enddo
    end do


    return
  end subroutine matvecp
!----------------------------
  subroutine vecvecp(current_mesh,vA,VB,Vprod)

    use global_data
    USE MYMPI

    implicit none

    TYPE(mesh_pointers) :: current_mesh
    real*8,POINTER ::  vA(:,:,:),VB(:,:,:)
    real*8,Intent(OUT) :: Vprod
!local variables
    Integer :: i,k,M,j,it
    real*8 :: myvprod
!---------------------------------------

    myvprod = zero
    do j = 0,current_mesh%ye
       do k = 0,current_mesh%ze
          do i = 0,current_mesh%xe
             myvprod = myvprod + VA(i,k,j)*VB(i,k,j)
          enddo
       enddo
    enddo

    call mpi_allreduce(myVprod,Vprod,1,MPI_DOUBLE_PRECISION,MPI_SUM,current_mesh%comm3d,ierr)
    return
  end subroutine vecvecp

END MODULE mg_solver
