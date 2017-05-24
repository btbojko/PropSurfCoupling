!**************************************************************************
MODULE data_types


  IMPLICIT NONE

!**************************************************************************

 INTEGER, PARAMETER ::  double =  8, unit_out      = 6, &
                        ngh = 2                         !number of ghost cells

 REAL *8, PARAMETER ::  zero    = 0.d0,   &
                        one     = 1.0d0,  &
                        two     = 2.0d0,  &
                        three   = 3.0d0,  &
                        four    = 4.0d0,  &
                        eight   = 8.0d0,  &

                        half    = 0.5d0,        &
                        quarter = 0.25d0,       &
                        eighth  = 0.125d0,      &
                        third   = one/three,         &
                        two_thirds = two/three,      &
                        root_two = 1.41421356237309515d0,   &
                        prec_ppe = 1.d-2,                   &
                        prec_gauss = 1.d-3,                   &

                        j2cal = 1.0d0/4.1868d0,         &
                        kg2g  = 1000.0d0,               &
                        m2cm  = 100.0d0,                &
                        runiv    = 1.9859d0,            &
                        mpa2atm = 9.86923266716d0,      &
                        pa2atm = mpa2atm*1d-6,          &     
                        kgmc2gcc = kg2g/m2cm**3,        &
                        j_msq2cal_cmsq = j2cal/m2cm**2, &
                        gcc2kgmc = one/kgmc2gcc,        &
                        kcalmc2atm = 4186.8/101325.d0,  &
                        j_kg2cal_g = j2cal/kg2g
                        
!!!SELECTED_REAL_KIND(P=14,R=30)  ! selects precision of real variables

!***************************************************************************************************

! Derived type containing global vectors used in the multigrid solution algorithm

  TYPE global_vectors

    REAL(KIND=double), DIMENSION(:), POINTER :: r
    REAL(KIND=double), DIMENSION(:), POINTER :: kp
    REAL(KIND=double), DIMENSION(:), POINTER :: p

    REAL(KIND=double), DIMENSION(:),  POINTER :: x1
    REAL(KIND=double), DIMENSION(:),  POINTER :: x2
    REAL(KIND=double), DIMENSION(:),  POINTER :: residual_vector
  
  END TYPE global_vectors

!***************************************************************************************************

! Derived type containing pointers for a given mesh

  type coarse_to_fine
     REAL(KIND=double) :: alpha(3)
     integer :: indx(3)
  end type coarse_to_fine
  type GMvec
     REAL(KIND=double), DIMENSION(:,:,:), POINTER :: mv
  end type GMvec
  type cppevec
     REAL(KIND=double), DIMENSION(:,:), POINTER :: ccmg
  end type cppevec
  TYPE mesh_pointers

    INTEGER :: mesh_num  ! number of mesh (1=finest,NOMSHS=coarsest)

    INTEGER :: xe,ze,ye  ! number of points in coordinate direction
    INTEGER :: nclr(2),nsnd(4,2),nrcv(4,2)
!
!   These are for the situation in which a processor has 0 points
!
    INTEGER :: myid,dnbr(4),comm3d,send_void(4),recv_void(4),rowcomm(2),rowid(2),rowsize(2)
    LOGICAL :: pskip(2),blank

!
!   First index is for the node info, second is the location of the block,
!   third is the color
!
    INTEGER, DIMENSION(:,:), POINTER :: iclr,kclr
    INTEGER, DIMENSION(:,:,:), POINTER :: sndclr,rcvclr

    REAL(KIND=double) :: dx,dz,dy,c1

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: f,ff       ! pointer to vector F for the mesh
    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: x,xold           ! pointer to vector X for the mesh

    REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: w          ! pointer to unknown  rank 4
    REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: r         ! pointer to residual rank 4
    REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: cm        ! pointer to coeffmat rank 4
    REAL(KIND=double), DIMENSION(:,:,:,:,:), POINTER :: dr        ! pointer to jacobian rank 4

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: cd
    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: c3
    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: c7
    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: cj

    type(cppevec), DIMENSION(:,:), POINTER :: CPPE

    REAL(KIND=double), DIMENSION(:), POINTER :: y,ya,eta
    REAL(KIND=double), DIMENSION(:), POINTER :: ey
    REAL(KIND=double), DIMENSION(:), POINTER :: eyy
    REAL(KIND=double), DIMENSION(:), POINTER :: eya
    REAL(KIND=double), DIMENSION(:,:), POINTER :: chi
    REAL(KIND=double), DIMENSION(:,:), POINTER :: chipr

    type (coarse_to_fine), DIMENSION(:,:,:,:), POINTER :: c2f

    type (GMVEC), DIMENSION(:), POINTER :: GM


    TYPE(mesh_pointers), POINTER :: coarse_mesh  ! pointer to coarse mesh
    TYPE(mesh_pointers), POINTER :: fine_mesh    ! pointer to fine mesh

    
  END TYPE mesh_pointers


  type  fourier_matrix
     integer :: nFourier
     real*8,POINTER :: A(:),M(:),F(:)
  end type  fourier_matrix
  
  type probe_type
     integer :: i,k,j
     real*8 :: T,X,Z,Y
     real*8 :: phi
  end type probe_type

!***************************************************************************************************


END MODULE data_types

!***************************************************************************************************
