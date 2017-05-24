MODULE GLOBAL_DATA

  USE data_types

  IMPLICIT NONE

  SAVE

  INTEGER, PARAMETER :: maxpack=3600, maxcells=100, ndim=3, maxproc = 500
  INTEGER, PARAMETER :: n_surf_terms = 4, num_rxn = 4, num_gas = 6, maxsize=7

  INTEGER :: neqmax, neqgas, irxn, iwflag
  INTEGER :: nn,nx,ny,ny_l_re
  INTEGER :: nyprocs,mywid,doparmdump,doprogress
  INTEGER :: npts,iord,ipack,ndump,frame,irestart,dooutput,doquaddump
  INTEGER :: nxproc,nzproc,nproc,ierr,nz,maxrange,flagfile,premixed
  INTEGER :: writemod,dflag,ncircles,drange(6),dnbr(6),xr(3),dotiming,myid
  INTEGER :: rowcom(2),colcom,comm3d,coords(3),pdims(3),surfid,surfcom,ddrange(6),rowsize(2)
  INTEGER :: rmap(4),mg_max4,n_proj, itime_sos
  INTEGER :: iperiod, nperiod, eqn_ref
  INTEGER :: ncyc,ncyc_init,ncyc_oseen,n_print_terms,ncyc_steady,ncyc_constant
  INTEGER :: bind_sizes,hb,mflux
  
  LOGICAL :: isperiodic(3),remain(3),isparallel(3),converged,issymmetric,isoddcoarse,ismgtxy
  LOGICAL :: doprojection,is2D, is_firstorder,do_solve_direct,output_bin

  REAL*8  :: tinit,tstop,time_steady_pressure,time_oseen
  REAL*8  :: alpha,beta,beta_pack,AP_frac,rho1
  REAL*8  :: xend,yend,dx,dy,c1y,dlength,timestep,pi,xloc,zend
  REAL*8  :: cp,pr,tcold,theta_ap,theta_binder,dz,lambda_eff
  REAL*8  :: da_ap,speed,da_binder,qheat_ap,qheat_binder
  REAL*8  :: rho_ap,rho_binder,rho_premixed
  REAL*8  :: lambda_ap,lambda_binder,pri_fac
  REAL*8  :: press,epi_p,omg_p,dpress_dt,phase_p

  REAL*8  :: storage_term,storage_term_2,storage_term_3
  REAL*8  :: period_p,cfl,period,dlx,dly,dlz,alp_W,alp_V,lewis(num_gas)
  REAL*8  :: tcyc,pcyc,coe_lamb1,coe_lamb2
  REAL*8  :: period_begin,phi_begin,  t_frame, phi_per, t_per, mass_flux
  REAL*8  :: starttime,commtime(2),gausstime(2), maxres(3),conv_error
  REAL*8  :: zplot, dt, dtmin, dt_timeint, dt_surface
  REAL*8  :: coe_d, coe_dt, time_coeff 
  REAL*8  :: alpha_test
  REAL*8  :: bind_diam(5),bind_mass(5)


  INTEGER,ALLOCATABLE :: nquad(:),quad(:,:),nxquad(:,:),nzquad(:,:)
 
  REAL*8,ALLOCATABLE :: tmdpnt(:,:,:),rhog(:,:,:),rhos(:,:,:),rhos_old(:,:,:)
  REAL*8,ALLOCATABLE :: qheats(:,:,:),lambdas(:,:,:),da_s(:,:,:)
  REAL*8,ALLOCATABLE :: theta_s(:,:,:),xsurf(:,:,:)
  REAL*8,ALLOCATABLE :: lambdag(:,:,:),dlgdx(:,:,:)
  REAL*8,ALLOCATABLE :: dlgdy(:,:,:),dlgdz(:,:,:),phi(:,:),dphidx(:,:),dphidz(:,:)
  REAL*8,POINTER,DIMENSION(:,:) :: phit,ff,tab,bc_ptr,dxm,dxp,dym,dyp,L_phi,xy
  REAL*8,ALLOCATABLE :: dphi2dx2(:,:),dphi2dz2(:,:),oldphi(:,:),psi(:,:,:)
  REAL*8,ALLOCATABLE :: rb(:,:),ft0(:,:),d_rb(:,:),d_ft0(:,:),d_vterm(:,:)
  REAL*8,ALLOCATABLE :: xl(:),xm(:),zl(:),zm(:),xcenter(:),ycenter(:),itype(:)
  REAL*8,ALLOCATABLE :: rad(:),x(:),y(:),z(:),detady(:),deta2dy2(:)
  REAL*8,ALLOCATABLE :: slambda(:),rlambda(:),sbuf1(:),sbuf2(:),rbuf1(:),rbuf2(:)
  REAL*8,ALLOCATABLE :: myright_ff_sbuffer(:),myright_ff_rbuffer(:)
  REAL*8,ALLOCATABLE :: myleft_ff_sbuffer(:),myleft_ff_rbuffer(:)
  REAL*8,ALLOCATABLE :: phi_movie(:),z_movie(:)

  CHARACTER(LEN=21) :: errfile,statfile,restartfile,outputfile

!***************************************************
!IMPLICIT VARIABLES (JACOBIANS)
!***************************************************
  INTEGER, PARAMETER :: max_new = 6, max_gss = 5;
  INTEGER :: it_new,it_gss
  REAL*8, PARAMETER ::  tol_gss = 1.d-6, omg=1.10d0
  REAL*8 :: res_gas, res_sld
  REAL*8,ALLOCATABLE,TARGET :: Bm(:,:,:,:,:)
  REAL*8,ALLOCATABLE :: oldrate(:,:,:),dtx(:,:,:)

!***************************************************
!FLUIDS VARIABLES
!***************************************************
  TYPE(mesh_pointers), POINTER :: finest_mesh    ! pointer to finest mesh data
  TYPE(mesh_pointers), POINTER :: coarsest_mesh  ! pointer to coarsest mesh data
  REAL(KIND=double), DIMENSION(:,:,:), POINTER :: wrk
  REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: wrk4


  REAL*8  :: dim_fact
!
  REAL*8,ALLOCATABLE :: ya(:),detadya(:)
  REAL*8,ALLOCATABLE :: vec(:,:),dphidxa(:,:),dphidza(:,:)
  REAL*8,ALLOCATABLE :: vvel(:,:,:),uvel(:,:,:),wvel(:,:,:),vcnv(:,:,:)
  REAL*8,ALLOCATABLE,TARGET :: p(:,:,:),divt(:,:,:),pold(:,:,:),wrk_vec(:,:,:,:)
  REAL*8,ALLOCATABLE,TARGET :: coeff(:,:,:,:),dfdt(:,:,:,:),rate(:,:,:,:),csld(:,:,:,:)
  REAL*8,POINTER :: q(:,:,:,:),f(:,:,:,:),oldsoln(:,:,:,:)
  REAL*8,ALLOCATABLE :: lbdavg(:,:,:,:),fdiff(:,:,:,:),dqdt(:,:,:,:) 
!************************************************************************
! MG PARAMETERS
!
!************************************************************************

  INTEGER :: gamma,num_mg_meshes,num_meshes,iguess,icoarse
  INTEGER :: mg_maxcycle,nu1,nu2,nu_fine,icycle,gauss_cycles

  INTEGER :: nnxx,nnyy,lvlx,lvly,ixp,jyq,iex,jey   !z and k are the same

  INTEGER :: mg_prm(4)

  INTEGER,ALLOCATABLE :: nx_mg(:),nz_mg(:),ny_mg(:),nx_tot(:),nz_tot(:),lx_mg(:)
  INTEGER,ALLOCATABLE :: ib_mg(:),kb_mg(:),id_mg(:),kd_mg(:)

  REAL(kind=double),ALLOCATABLE :: dx_mg(:),dz_mg(:),dy_mg(:)

  REAL(kind=double),ALLOCATABLE :: dphidx_mg(:,:,:),dphidz_mg(:,:,:)
  REAL(kind=double),ALLOCATABLE :: dphidxa_mg(:,:,:),dphidza_mg(:,:,:)

  LOGICAL,ALLOCATABLE :: iscoarse_xz(:), iscoarse_y(:)


  INTEGER  :: nreaction, nclip
  REAL*8   :: qg(num_rxn)
  REAL*8   :: theta2_2_theta3,pre2_2_pre3
  REAL*8   :: rho_m_st,cp_ap,cp_binder
  REAL*8   :: mu(num_rxn,num_gas),yf(num_gas,num_rxn)
  REAL*8   :: alpha_pack, alpha_fdf,th_st, th_ap, th_bind,th_final
  REAL*8   :: Q_ap,Q_binder,A_ap,A_binder,E_ap,E_binder
  REAL*8   :: c_lcp,f_pyr,A_ox,E_ox
  REAL*8   :: surfvel_max, surfvel_min

  REAL*8,ALLOCATABLE :: ay(:),by(:),cy(:),zy(:),fy(:),fyv(:,:),Binv(:,:),Bb(:,:,:)
  REAL*8,ALLOCATABLE :: aq(:),cq(:),Vdiag0(:),Vright1(:)
  REAL*8,ALLOCATABLE :: thetay(:),da(:),np(:)


  REAL*8, POINTER  :: myphi(:),allphi(:),G_phi(:)
  REAL*8, POINTER  :: G_f(:,:,:),G_x(:),G_speed(:),G_rate(:,:,:),G_q(:,:,:)
  REAL*8, POINTER  :: G_psi(:,:)
!!
!!
  REAL*8   ::  rb_melt, T_melt, T_bmelt, g_accel(ndim)
  INTEGER  ::  i_flame_pos

! SAVED VARIABLES IN LU DECMPOSITION LOOPS FOR TXY (this is not the sparse LU solver)
! these should save you time but increase memory requirement
  REAL*8,ALLOCATABLE :: L_save(:,:,:,:,:),U_save(:,:,:,:,:),Binv_save(:,:,:,:,:)
  REAL*8,ALLOCATABLE ::  Lss(:,:,:),Uss(:,:,:)
  LOGICAL :: saved_LU_vars
  REAL*8 :: times(2,20)

!VARIABLES INVOLVED IN THE BOUNDARY CONDITIONS
  LOGICAL :: is_proc_bnd(4),is_periodic_xz(0:4),is_periodic_y, use_extrap(0:3,4)
  LOGICAL :: time_varying_inflow, is_TVD, f1st_ord_conv, is_sandwich
  integer :: nxv(0:ndim+1),nyv(0:ndim+1),ibc(0:4,4,2),sz_bndry
  integer :: eqn_nclr(2,0:3),addptx(0:4),addpty(0:4)
  INTEGER, ALLOCATABLE :: type_bc(:,:,:)
  INTEGER,ALLOCATABLE :: eqn_iclr(:,:,:),eqn_kclr(:,:,:)
  REAL*8,ALLOCATABLE :: vec_dirichlet(:,:,:), v_south(:),qsurf(:,:,:)
  INTEGER,ALLOCATABLE :: nx_mg_fix(:),nx_tot_fix(:)
  REAL*8   :: addbnd(1:4), xstart_MG, Q_scalar(ndim,4)

!VARIABLES for the new (Tom's) packs
  REAL*8   :: aspect_ratio,phi0
  INTEGER  :: iper
!VARIABLES for avg
  REAL*8   :: surf_terms(n_surf_terms), surf_terms0(n_surf_terms)
  REAL*8   :: rb0save,rb00,dt_integ,time_period,pos_period
  REAL*8   :: rblend,lblend,alfa_blend
  INTEGER  :: ndump_s_terms

!some new stuff
  REAL*8,ALLOCATABLE  :: dconv(:,:,:)
  INTEGER,ALLOCATABLE :: iconv(:,:,:)
  LOGICAL :: converged_split(maxsize)
  INTEGER  :: var_updated
  CHARACTER*28 :: flag_comp

END MODULE GLOBAL_DATA
!********************************************************************************

!**************************************************************
MODULE LU_VARS

  SAVE
!********************************
!   SPARSE LU
!********************************
  integer :: LU_nx, LU_ny, nmaxj, nmaxi, ndiag, sizea, LU_frequency
  real*8,ALLOCATABLE  :: as(:,:),ab(:,:), ls(:,:),ub(:,:),us(:,:)
  real*8,ALLOCATABLE  :: LU_y(:),LU_f(:)
  real*8,ALLOCATABLE  :: cm(:,:,:)
  real*8,ALLOCATABLE  :: buff(:),buff2(:),buff3(:)
  real*8              :: residual_direct
  integer,ALLOCATABLE :: ja(:,:),ia(:,:),nja(:),nia(:)
  integer,ALLOCATABLE :: mainD_j(:),mainD_i(:)
  integer,ALLOCATABLE :: njl(:),jl(:,:),niu(:),iu(:,:)
  integer,ALLOCATABLE :: nju(:),ju(:,:)
  integer,ALLOCATABLE :: posdiag(:,:)
  integer,ALLOCATABLE :: LU_size(:), LU_disp(:)
  integer,ALLOCATABLE :: allcounts(:), displacements(:)
  integer :: lastmem, LU_last
  integer :: LU_prd(2),LU_np(2), LU_diag_gather
  logical :: done_LU_decomposition

END MODULE LU_VARS
!********************************************************************************
