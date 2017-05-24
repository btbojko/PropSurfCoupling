!     *********************************************************************************
!     **********************************MULTIGRID**************************************

MODULE mg_solver

  USE data_types
  USE global_data
  USE LUSOLVER

  IMPLICIT NONE

  PUBLIC  :: fmg_solve           ! solves the problem using the full multigrid method

  PRIVATE :: mg_cycle,        &  ! performs one complete multigrid cycle
             t_x_calc,        &  ! interpolates the correction from the coarse mesh to the fine mesh
             t_trans_r_calc      ! restricts the residual from the fine mesh to the coarse mesh

CONTAINS

!***************************************************************************************************

  SUBROUTINE fmg_solve

  ! Subroutine to perform complete full multigrid cycles on the matrix equation Kx = f
  ! until convergence is achieved.  The process is described in the flowchart on page 45
  ! of: "Multigrid Methods: Fundamental Algorithms, Model Problem Analysis and 
  ! Applications" by K.STUBEN and U.TROTTENBERG, in "Lecture Notes in Mathematics 960:
  ! Multigrid  Methods", edited by W.Hackbusch and U.Trottenberg, Springer-Verlag 1982.
  !
  ! Note: the finest mesh is mesh #1, 
  !       the coarsest mesh is mesh #NOMSHS.

!---------------------------------------------------------------------------------
!   Local variables:

    TYPE(mesh_pointers), POINTER :: current_mesh

    REAL(KIND=double) :: norm_f, norm_r, time

    INTEGER :: nxf,nzf,nyf,iter
!---------------------------------------------------------------------------------

    icycle = 0
!
!   Compute the load vector on each coarse mesh by restricting the load
!   

    current_mesh => finest_mesh

    IF (iguess.eq.1 .OR. num_mg_meshes == 1) THEN

       do iter = 1, mg_maxcycle
       CALL gauss_relax( current_mesh, nu_fine, .TRUE.)   !> USE .F. <
       enddo
       RETURN

    ELSEIF (iguess == 2) THEN

       converged = .FALSE.
       CALL gauss_relax( current_mesh, 2*nu_fine, .FALSE.)   !> USE .F. <

    ELSE


      DO

      IF ( .NOT. ASSOCIATED(current_mesh%coarse_mesh) ) EXIT

      nxf=current_mesh%xe +ngh;nzf=current_mesh%ze +ngh;nyf=current_mesh%ye +ngh
      wrk(-ngh:nxf,-ngh:nzf,-ngh:nyf) = current_mesh%f(-ngh:nxf,-ngh:nzf,-ngh:nyf)

!      CALL restrict( current_mesh%coarse_mesh )
      CALL transfer( current_mesh%coarse_mesh )

      current_mesh => current_mesh%coarse_mesh

      END DO

      IF(do_solve_direct) then
         CALL SOLVE_DIRECT( coarsest_mesh)
      ELSE
         CALL gauss_relax( coarsest_mesh, nu1+nu2, .FALSE. )  !USE > .F. <
      ENDIF

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

   END DO

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

        CALL gauss_relax( current_mesh, nu1, .TRUE. )   !>USE .T.<
 
        c(current_mesh%mesh_num) = c(current_mesh%mesh_num) + 1

        CALL restrict( current_mesh%coarse_mesh )

        current_mesh => current_mesh%coarse_mesh

        current_mesh%x = zero

        END DO
 

      IF(do_solve_direct) then
         CALL SOLVE_DIRECT( coarsest_mesh)
      ELSE
         CALL gauss_relax( coarsest_mesh, nu1+nu2, .FALSE. )  !USE > .F. <
      ENDIF

      current_mesh => current_mesh%fine_mesh
    
      IF ( .NOT. ASSOCIATED(current_mesh) ) RETURN

        DO

        CALL prolong( current_mesh%coarse_mesh, current_mesh )

        CALL gauss_relax( current_mesh, nu2, .FALSE. )    !>USE .F.<

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

!***************************************************************************************************

END MODULE mg_solver

!***************************************************************************************************
