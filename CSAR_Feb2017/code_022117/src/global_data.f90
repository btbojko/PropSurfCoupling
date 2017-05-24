
!********************************************************************************
MODULE MQDATA
  USE data_types

  IMPLICIT NONE

  SAVE

  REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:)  :: MQvelocity
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:)  :: MQterm
  REAL*8,ALLOCATABLE,DIMENSION(:,:)  :: MQeta,MQchi,MQchipr,MQchisc,MQphi
  REAL(kind=double),ALLOCATABLE,DIMENSION(:,:,:) :: phi_mg,psi_mg,temp_mg
  REAL*8 :: ystart,detabar,deta,MQgy,fmax,velofmax,MQnonC(8)
  REAL*8  :: dyterm(2),convTerm(3,2),uwvTerm(3),diffTerm(3,2)
  INTEGER :: quartuple(4),qpick(4)
  INTEGER :: quintuple(5),qpick5(5)
  INTEGER :: jstart,jend,mcomp,mgradj,mgradjM,Mindx

END MODULE MQDATA
!*********************************************************************************

MODULE GLOBAL_DATA

  USE data_types
  USE MQDATA

  IMPLICIT NONE

  SAVE

  INTEGER, PARAMETER :: maxsize=5
  INTEGER, PARAMETER :: ndim=3

  INTEGER, PARAMETER :: maxsizeGMRES = 23
  INTEGER, PARAMETER :: maxpack=5000, maxcells=10


  INTEGER :: nx,ny,neqmax,nyprocs,doparmdump,doprogress
  INTEGER :: npts,iord,iperiod,nperiod,ipack,ndump,nsfxrstrt,frame,frameCUT
  INTEGER :: irestart,dooutput,doquaddump
  INTEGER :: ierr,nz,method,maxrange,flagfile,ncuts,icut
  INTEGER :: writemod,dflag,ncircles,ncirclesIN
  INTEGER :: drange(6),dnbr(6),ncells,xr(3),dotiming,dwiderange(6)
  INTEGER :: comm3d,cutcomm,surfcom,cutxcomm,cutzcomm
  INTEGER :: nxproc,nzproc,nproc, nproc_cut, mywid,myid,cutid,cutXid,cutZid
  INTEGER :: coords(3),pdims(3),surfid,dSerialRange(6),rowsize(2)
  INTEGER :: rmap(4),mg_max4,kplot, itime_sos, done_period
  INTEGER :: jRFLO(3)=0,maxFvars,maxQvars
  INTEGER :: ncyc,ncyc_init,ncyc_oseen, ncyc_run, n_samples,ncyc_debug
  INTEGER :: Tjumps,Tncyc,Tfrac
  INTEGER :: inewton,Nnewton
  INTEGER :: ntrck
  INTEGER :: nsteps
  INTEGER,allocatable :: itrck(:),jtrck(:),ktrck(:),mtrck(:)   !locals for now


  LOGICAL :: isperiodic(3),remain(3),isparallel(3),converged,issymmetric,isoddcoarse,ismgtxy
  LOGICAL :: doprojection, is_firstorder, is2D, converged_split(maxsize), skipTVD,startVIZ2
  LOGICAL :: skipDeformation,skipTemperature,skipViscDeriv,skipContSource,skipPressure
  LOGICAL :: skipDeformationFluid,SkipVarDensity,skipCrossFlow,skipUnsteady,skipNorthNeumanBcTemp
  LOGICAL :: SameTopology,steadyconvergence = .false.,numerdivergence  = .false.

  REAL*8  :: tinit,xend,yend,dx,dy,c1y,dlength,tstop,timestep,time_cyc,timeStepOld
  REAL*8  :: alpha,cp,pr,equivalence_ratio,beta,pi,xloc,zend
  REAL*8  :: tcold,theta_ap,theta_binder,dz,lambda_eff,rho_eff
  REAL*8  :: da_ap,da_binder,qheat_ox,qheat_binder,rho_ap,rho_binder
  REAL*8  :: lambda_ap,lambda_binder,press,epi_p,omg_p
  REAL*8  :: period_p,cfl,period,periodXZ,dlx,dly,dlz,lewis(maxsize),tcyc(2)
  REAL*8  :: alp_W(maxsize),alp_V(maxsize)
  REAL*8  :: period_begin,phi_begin,  t_frame, phi_per, t_per, mass_flux
  REAL*8  :: starttime,commtime(2),gausstime(2), elapsedTime,maxtime,done_time
  REAL*8  :: zplot,xplot, dt, dtmin,Tlimit(2),maxres(4),conv_error
  REAL*8  :: coe_d, coe_dt, time_coeff,coe_exp
  REAL*8  :: T_melt=834d0, rb_melt, TRFLO,rhoRef
  REAL*8  :: speed,speed_timeint=0d0, dt_timeint
  REAL*8  :: rho1,alpha_2dpack
  REAL*8  :: da(maxsize),thetay(maxsize),zetay(maxsize),np(maxsize),qg(maxsize)
  REAL*8  :: lengthscale

  !>>> vic: shear, from juPack
  REAL*8  :: shear_rate, shearFreq, shearRateAmp
  !<<< vic

  INTEGER,ALLOCATABLE :: nquad(:),quad(:,:),nxquad(:,:),nzquad(:,:)
 
  REAL*8,ALLOCATABLE :: tmdpnt(:,:,:),rhos(:,:,:)
  REAL*8,ALLOCATABLE :: qheats(:,:,:),lambdas(:,:,:),da_s(:,:,:)
  REAL*8,ALLOCATABLE :: theta_s(:,:,:),xsurf(:,:,:),psurf(:,:,:)
  REAL*8,ALLOCATABLE :: lambdag(:,:,:),rhog(:,:,:)
  REAL*8,POINTER     :: phi(:,:)
  REAL*8,ALLOCATABLE :: dlgdx(:,:,:),dlgdy(:,:,:),dlgdz(:,:,:)
  REAL*8,ALLOCATABLE :: dphidx(:,:),dphidz(:,:),phit(:,:)
  LOGICAL,ALLOCATABLE :: taGsmooth(:,:)
  REAL*8,ALLOCATABLE :: dphi2dx2(:,:),dphi2dz2(:,:),oldphi(:,:),psi(:,:,:)
  REAL*8,ALLOCATABLE :: vec(:,:,:),dphidxa(:,:),dphidza(:,:),vec0(:,:)
  REAL*8,ALLOCATABLE :: rb(:,:),ft0(:,:),d_rb(:,:),d_ft0(:,:),d_vterm(:,:)
  REAL*8,ALLOCATABLE :: xl(:),xm(:),zl(:),zm(:),xcenter(:),ycenter(:),zcenter(:)
  REAL*8,ALLOCATABLE :: rad(:),x(:),ya(:),z(:),detady(:),deta2dy2(:)
  REAL*8,ALLOCATABLE :: slambda(:),rlambda(:),sbuf1(:),sbuf2(:),rbuf1(:),rbuf2(:)
  REAL*8,ALLOCATABLE :: myright_ff_sbuffer(:),myright_ff_rbuffer(:)
  REAL*8,ALLOCATABLE :: myleft_ff_sbuffer(:),myleft_ff_rbuffer(:)
  REAL*8,ALLOCATABLE :: phi_movie(:),z_movie(:)
  REAL*8, ALLOCATABLE :: dconv(:,:,:,:)
  INTEGER,ALLOCATABLE :: iconv(:,:,:,:)
  Character(5),ALLOCATABLE :: OxidizerType(:)

!these are time extrapolated values
  !REAL*8,POINTER :: Edphidx(:,:),Edphidz(:,:),Edphi2dx2(:,:),Edphi2dz2(:,:)
  !REAL*8,POINTER :: Edphidxa(:,:),Edphidza(:,:)
  REAL*8,POINTER :: y(:)

  CHARACTER(LEN=21) :: errfile

!
! Filename variables
!
  CHARACTER :: inputfile*100
  CHARACTER :: inputfile_name*100
  CHARACTER :: chem_inputfile_name*100
  CHARACTER :: pack_inputfile_name*100
  INTEGER :: inreadfilename

!***************************************************
!IMPLICIT VARIABLES (JACOBIANS)
!***************************************************
  INTEGER, PARAMETER :: max_new = 6, max_gss = 5;
  INTEGER :: it_new,it_gss
  integer :: gmrssize,gmrssize2
  REAL*8, PARAMETER ::  tol_gss = 1.d-6, omg=1.10d0
  REAL*8 :: res_gas, res_sld,errGMres
  REAL*8,ALLOCATABLE,TARGET :: Bm(:,:,:,:,:),vgm(:,:,:,:,:)
  REAL*8,ALLOCATABLE :: dgm(:),hgm(:,:),c1gm(:),c2gm(:),ygm(:)
  REAL*8,ALLOCATABLE :: oldrate(:,:,:),dtx(:,:,:),rhsold(:,:,:,:),dtorhocp(:,:,:)

!***************************************************
!FLUIDS VARIABLES
!***************************************************
  TYPE(mesh_pointers), POINTER :: finest_mesh    ! pointer to finest mesh data
  TYPE(mesh_pointers), POINTER :: coarsest_mesh  ! pointer to coarsest mesh data
  REAL(KIND=double), DIMENSION(:,:,:), POINTER :: wrk,pold,p
  REAL(KIND=double), DIMENSION(:,:,:,:), POINTER :: wrk4


  REAL*8  :: dim_fact,factorRB
!
  REAL*8,ALLOCATABLE :: detadya(:),deta2dy2a(:)
  REAL*8,ALLOCATABLE :: vvel(:,:,:),uvel(:,:,:),wvel(:,:,:),vcnv(:,:,:),d_vvel(:,:,:)
  REAL*8,POINTER :: divt(:,:,:),Postdivt(:,:,:)
  REAL*8,ALLOCATABLE,TARGET :: coeff(:,:,:,:,:),dfdt(:,:,:,:),rate(:,:,:,:),csld(:,:,:,:),inirate(:,:,:,:)
  REAL*8,ALLOCATABLE :: ccs1(:,:,:),ccs2(:,:,:),ccs3(:,:,:),ccs4(:,:,:),ccs5(:,:,:),ccs6(:,:,:)
  REAL*8,ALLOCATABLE :: ccs7(:,:,:),ccs8(:,:,:),ccs9(:,:,:)
  REAL*8,POINTER :: q(:,:,:,:),f(:,:,:,:),oldsoln(:,:,:,:),inif(:,:,:,:)
  REAL*8,ALLOCATABLE :: newQBnd(:,:,:),srfQBnd(:,:,:)
  REAL*8,ALLOCATABLE,TARGET :: dff(:,:,:,:)
  REAL*8,ALLOCATABLE :: lbdavg(:,:,:,:),dqdt(:,:,:,:)  !maybe you can drop it
  REAL*8,ALLOCATABLE :: rate_out(:,:,:,:)

  !>>> vic: shear, from juPack
  REAL*8, ALLOCATABLE :: shear(:)
  !<<< vic

!************************************************************************
! MG PARAMETERS
!
!************************************************************************

  INTEGER :: gamma,num_mg_meshes,num_meshes,iguess,icoarse
  INTEGER :: mg_maxcycle,nu1,nu2,nu_fine,icycle

  INTEGER :: nnxx,nnyy,lvlx,lvly,ixp,jyq,iex,jey   !z and k are the same

  INTEGER :: mg_prm(4)

  INTEGER,ALLOCATABLE :: nx_mg(:),nz_mg(:),ny_mg(:),nx_tot(:),nz_tot(:),lx_mg(:)
  INTEGER,ALLOCATABLE :: ib_mg(:),kb_mg(:),id_mg(:),kd_mg(:),gmres_size_mg(:)

  REAL(kind=double),ALLOCATABLE :: dx_mg(:),dz_mg(:),dy_mg(:)

  REAL(kind=double),ALLOCATABLE :: dphidx_mg(:,:,:),dphidz_mg(:,:,:)
  REAL(kind=double),ALLOCATABLE :: dphidxa_mg(:,:,:),dphidza_mg(:,:,:)

  LOGICAL,ALLOCATABLE :: iscoarse_xz(:), iscoarse_y(:)


  INTEGER  :: neqgas,irxn,nspecies
  INTEGER  :: nclip(maxsize+2)=0
  INTEGER  :: nclipwrite
  REAL*8   :: n1,n2,n3,n4, theta1,theta2,theta3,theta4
  REAL*8   :: da1,da2,da3,da4, qg1,qg2,qg3,qg4
  Real*8   :: z1,z2,z3,z4
  REAL*8   :: theta2_2_theta3,pre2_2_pre3
  REAL*8   :: rho_m_st,cp_ap,cp_binder,cp_gas
  REAL*8   :: mu(maxsize+2,maxsize+2),yf(maxsize+2,maxsize+2)
  REAL*8   :: rhoV_RocPack, alphaW_Rocpack,beta_pack
  REAL*8   :: alpha_pack,rho_pack,alphaH_pack,rhoH_pack,lambda_pack
  REAL*8   :: th_st, th_ap
  REAL*8   :: Tflame1, Tflame2,TheoreticalPackingDensity
  REAL*8   :: Q_ap,Q_binder,A_ap,A_binder,E_ap,E_binder
  REAL*8   :: A_ox,E_ox
  REAL*8   :: surfvel_max, surfvel_min
  REAL*8   :: a_lamg,b_lamg,e_lamg,c_lcp,exp_lcp
  REAL*8   :: a_cpg,b_cpg,e_cpg,T_cpg

  REAL*8,ALLOCATABLE :: ay(:),by(:),cy(:),zy(:),fy(:),fyv(:,:),Binv(:,:),Bb(:,:,:)
  REAL*8,ALLOCATABLE :: aq(:,:),cq(:,:),Vdiag0(:),Vright1(:)


!visualization // ALLGATHER variables
  REAL*8, POINTER  :: myphi(:),allphi(:),G_phi(:,:)
  REAL*8, POINTER  :: G_f(:,:,:),G_x(:),G_speed(:),G_rate(:,:,:),G_q(:,:,:,:)
  REAL*8, POINTER  :: G_psi(:,:),G_z(:)
  

! SAVED VARIABLES IN LU DECMPOSITION LOOPS FOR TXY (this is not the sparse LU solver)
! these should save you time but increase memory requirement
  REAL*8,ALLOCATABLE :: L_save(:,:,:,:,:),Binv_save(:,:,:,:,:)
  LOGICAL :: saved_LU_vars

!LOCATIONS
  REAL*8,POINTER :: yloc(:), myphi1(:),allphi1(:),myphi2(:),allphi2(:),myphi3(:),allphi3(:)
  INTEGER,ALLOCATABLE :: jloc(:)
  INTEGER, ALLOCATABLE  :: allsize(:),displacements(:),alldrange(:),alldrange11(:,:)
  INTEGER ,ALLOCATABLE :: allsize1(:), displacements1(:),alldrange1(:)
  INTEGER ,ALLOCATABLE :: allsize2(:), displacements2(:),alldrange2(:)
  INTEGER  ::  nyloc,inyloc

!ADDED

  INTEGER       :: nxv(0:5),nzv(0:5),nyv(0:5),n_print_terms,ncyc_steady
  INTEGER       :: iwflag,mysizeB,mysizeY
  CHARACTER*25  :: restartfile,outputfile,flag_comp
  REAL*8        :: pcyc, storage_term_2,phase_p,time_steady_pressure,time_oseen,g_accel(3)
  INTEGER       :: rowcomm(2),colcomm,rowid(2),colid
  character*3   :: sfx
  character*4   :: sfxrestart,sfxtop
  character*20  :: TimeString,TimeStringRestart

!added to module TRYEROS
  INTEGER      :: Nimposed
  REAL*8 :: MImposed
  REAL*8,ALLOCATABLE :: Uimposed(:),Zimposed(:)

!source term

 REAL*8,ALLOCATABLE  :: rhogasCont(:,:,:,:),dtCont(:)


 type(fourier_matrix) :: matrix
 
 integer :: Nprobes
 type(probe_type),allocatable :: solidProbes(:),OldSolidProbes(:)


  Real*8   :: nus_ap,nus_binder,f_pyr(4),df_pyr
  Real*8   :: QsolidPar(3),theta_decomp,Da_decomp,ndecomp,T_AP,P_AP
  Real*8   :: BinderThickness

  real*8 :: ts_alum,ts_AP,ts_BD
  real*8 :: ts_alum_H,ts_AP_H,ts_BD_H    !volume_{*,H}/volume_matrix
  real*8 :: ts_alum_HT,ts_AP_HT,ts_BD_HT !volume_{*,H}/volume_tot
  real*8 :: ws_alum_H,ws_AP_H,ws_BD_H

  LOGICAL :: ICDEBUG

  INTEGER :: iloc,ploc,ncyc_p,NCYCRAD

  REAL*8,POINTER     :: Gfield(:,:,:),radheat(:,:,:),BCradheat(:,:)
  REAL*8,POINTER     :: VFrac(:,:,:,:),DiamP(:,:,:,:),Vod(:,:,:,:)
  real*8 :: beta_al,rho_Al,lambda_al,rho_Ala,lambda_ala,Initial_diameter,beta_rhoal
  real*8 :: max_f(maxsize),min_f(maxsize),KrierFACTOR

  real*8,pointer :: x_mg(:,:),z_mg(:,:)
  real*8,pointer :: deltaKron(:,:)
  character(LEN=10) :: sCASE
  logical :: NORADIATION,PRESSVARY,HEATVARY,NOQ4
  real*8 :: InitialPressure,initialTimestep

END MODULE GLOBAL_DATA
