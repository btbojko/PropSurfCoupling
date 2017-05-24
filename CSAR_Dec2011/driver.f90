!    *****************************************************************
!                          ROCFIRE3D START POINT 
!    *****************************************************************
PROGRAM ROCFIRE  
USE GLOBAL_DATA
USE BC_DATA
USE MYMPI
  
IMPLICIT NONE

include "parallel.h"

!--------------------------------------------------------------------
! Local variables
INTEGER ::  filenum, i, i5, m, n,j,k,idump,itemp
INTEGER ::  ii,imax,eqn,n_iter_surf
INTEGER ::  sum1,sum2,sum_K
REAL*8  ::  tout, t
REAL*8  ::  error, u3, v3, w3
REAL*8  ::  vol,vol_tot,vol_ap,vol_bind,vol_pack
REAL*8  ::  lengthscale, ratio
REAL*8  ::  ddlx,ddly,ddlz
REAL*8  ::  rvtot_K,rvexp_K,beta_K,diff_beta_K
REAL*8  ::  th,lb1,lb3,coe_rho
REAL*8  ::  t1,s1,s2,s3,t2,gm_ap,gm_tot
!
! PRF: Timing variables
!
REAL*8  ::  endtime,inputtime,ictime,quad2time
REAL*8  ::  packtime,setuptime,gridtime,udpsitime,iudpsitime
REAL*8  ::  spherebctime,quadtime,calctime,optime
REAL*8  ::  txytime, phitime,bctime,restarttime
REAL*8  ::  progtime,temptime,psitime,ictime1,maxsetuptime
REAL*8  ::  outtime,decomptime,total_timestep,maxtimestep
REAL*8  ::  PRESSURE

LOGICAL :: done_time

! PRF: Process specific output files
CHARACTER(LEN=21) :: rstartfile,headerfile,timingfile
CHARACTER(LEN=21) :: moviefile
!--------------------------------------------------------------------------
  
! PRF: Required for all MPI apps
call MPI_INIT(ierr)
  
call MPI_COMM_RANK(MPI_COMM_WORLD,mywid,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

starttime = MPI_WTIME()

CALL INPUTS
ncyc = 0
tcyc = 0
pcyc = press
itemp = 6
hb = 22
mflux=23

!Set Sizes
neqgas = num_gas
neqmax = neqgas+1
irxn = num_rxn

CALL ALLOCATE_INITIAL
  
IF (ipack <= 0) THEN        !Steady-state or sandwich calculation
    dlength = xend
    period = dlength
    dly = dlength

    alp_W = alpha       !mass fraction AP in binder
    alp_V = alp_W/(rho_ap/rho_binder*(1.0d0-alp_W)+alp_W)   !volume fraction AP in binder
  
    if(premixed == 0)then
        rho_premixed = rho_binder
    else
        rho_premixed = rho_ap*alp_V + rho_binder*(1.0d0-alp_V)
    endif
    rho1 = zero

    do_solve_direct = .TRUE.
        
    alpha_pack = alp_W
    if(ipack == 0) then  !!!>>sandwich<<!!!
        output_bin = .TRUE.
        is_sandwich = .TRUE.
        issymmetric = .TRUE.
        coe_d = one     !unsteady solid
        xloc = AP_frac*dlength
        gm_ap = AP_frac*rho_ap+(1.0d0-AP_frac)*rho_premixed*alp_W
        gm_tot = AP_frac*rho_ap+(1.0d0-AP_frac)*rho_premixed
        alpha_pack = gm_AP/gm_tot
        alpha_fdf = AP_frac*rho_ap/((1.0d0-AP_frac)*rho_premixed + AP_frac*rho_ap)
!        call read_chem_info
    elseif (ipack < 0) then    !steady 1-d computations
        xloc = 1.d10
        coe_d = zero            !steady solid
        ncyc_oseen = npts+2
        epi_p = zero
        nnxx = 1
        lvlx = 1
        writemod = 100
        irestart = 0
        alpha_pack = alp_W
        alpha_fdf = alp_W
        if(nproc >1)then
           write(*,*)'CASE ipack = -1 runs in single proc mode only'
           stop 'driver.f90'
        endif
    endif
ELSE
    output_bin = .TRUE.
    issymmetric = .FALSE.    
    do_solve_direct = .FALSE.
    alp_W = alpha       !mass fraction AP
    alp_V = alp_W/(rho_ap/rho_binder*(1.0d0-alp_W)+alp_W)   !volume fraction AP
    if(premixed == 0)then
        rho_premixed = rho_binder
    else
        rho_premixed = rho_ap*alp_V + rho_binder*(1.0d0-alp_V)
    endif
    CALL READ_DISCS
    coe_d = one
    xend = dlength
    zend = dlength 
ENDIF

! HARDWIRED FLAG ---------------------------
!
is_firstorder = .TRUE.
is2D = .TRUE.
!
is_periodic_xz = .TRUE.
is_periodic_y = .FALSE.
if(is_sandwich) then        
    is_periodic_xz = .TRUE.
    is_periodic_y = .FALSE.
endif
!
time_varying_inflow = .FALSE.
is_TVD =  .FALSE.
f1st_ord_conv = .FALSE.  !first order CONVECTION if not using a TVD scheme turn this flag on when the cross-flow is high
  
CALL run_parameter
CALL global_parameter
!
! ASSIGN number of processor in each dimension: pdims(1:3)    
!
CALL GET_PROC_DECOMP
  
CALL GET_MG_DECOMP  
  
CALL ALLOCATE_MAIN 

! Assign ouput file names
flagfile=press

WRITE(statfile,'("stats",I3.3,".dat")') flagfile
WRITE(errfile,'("history",I3.3,".dat")') flagfile
WRITE(timingfile,'("timing",I3.3,".dat")') flagfile
if(myid==0) then 
    open(Unit=hb,FILE='BurningRate.csv',status='unknown')
    write(hb,*) ' Cycle,       dt(s),     time(s),  Average_y(cm),   Rate(cm/s), Ave_Rate(cm/s)'  
    open(unit=mflux,File='Flux_avg.csv',status='unknown')
    write(mflux,*)'Cycle,     time(s),  Ave_mass_flux,   Mass_Flux2,mgas'
endif

!---- setup Grid and Initial Guess
CALL GRIDSETUP

!---- setup bc values, some may be dependent on the grid
CALL SET_CONDITIONS_BC  
  
IF (irestart == 0) THEN
    CALL INIT_FROM_SCRATCH
    if(myid == 0) open(UNIT=9,FILE=errfile,STATUS='UNKNOWN',POSITION='REWIND')
    if(myid == 0) open(UNIT=29,FILE=statfile,STATUS='UNKNOWN',POSITION='REWIND')
ELSEIF (irestart == 1) THEN
    CALL RESTART(restartfile)
    if(myid == 0) open(UNIT=9,FILE=errfile,STATUS='UNKNOWN',POSITION='REWIND')
    if(myid == 0) open(UNIT=29,FILE=statfile,STATUS='UNKNOWN',POSITION='REWIND')
ELSEIF (irestart == -1) THEN
    CALL RESTART_1D
    if(myid == 0) open(UNIT=9,FILE=errfile,STATUS='UNKNOWN',POSITION='REWIND')
    if(myid == 0) open(UNIT=29,FILE=statfile,STATUS='UNKNOWN',POSITION='REWIND')
ENDIF

!zero out non-state vars
uvel = zero;vvel =zero;wvel = zero
dtx = zero;coeff = zero;ft0 = zero;
    

!*******************************************************************************************
if (ipack == -1) then
    write(*,*)'***************************************'
    write(*,*)'----------------------------------------'
    write(*,*)
    write(*,*)'RUNNING with flag ipack= ',ipack
    write(*,*)'1D simulation<--->steady solid temperature'
    write(*,*)
    write(*,*)'----------------------------------------'
    write(*,*)'***************************************'
endif
!*******************************************************************************************
!      START COMPUTATION                     
  
times = zero;gausstime=zero;commtime=zero;
CALL PRN_GLOB_VARS

t = tinit
tout = t
tcyc = t

do nn=1,npts
    dt = timestep
    ncyc=ncyc+1
    tcyc = tcyc + timestep
    pcyc = pressure(tcyc)  !spatially constant
    CALL GET_TIME_COEFF 

!--- Update equations
    dt = timestep  ! if it changed during sub_cycling in TXY
    if( ipack == 1 ) CALL BC
    CALL UPDATE_TXY

    CALL CURB_TXY

    CALL DIRICHLET_BC(itemp,finest_mesh) 

!--- Update Boundary Conditions
    CALL BC
    CALL PYROLYSIS
    CALL STORAGE
   

!--- Update the explicit part of NS (Crank-Nicholson)
    IF (doprojection )THEN
        CALL BCQ
        CALL PROJECTION_SOURCE  !set rate to zero and rate(:,:,:,4)
        CALL EXP_RHSQ           !set rate(:,:,:,1:3)
    ENDIF
     
!--- Update surface
    n_iter_surf = 1
    if(ipack == 0 .and. ncyc < ncyc_steady) then
        CALL SANDWICH_PHI(t,tout)
        done_time = .FALSE.
    else
        CALL PROPAGATE_SURF(t,tout,n_iter_surf,done_time)
    endif

!
    if(doprojection) then
        CALL BCQ                   !Dirichlet Boundaries
        CALL UPDATE_UVW
        CALL PROJECTION
        CALL BCQ 
    endif        
    
    if(ncyc == 1309) then
        continue
    endif
     
    CALL PROGRESS(t,tout)

    t = tout
   
    CALL VISUALIZE_FIELD  

    if (ipack == 1) then
        iperiod = int(phi(0,0)/period) - int(phi_begin/period)
        call MPI_BCAST(iperiod,1,MPI_INTEGER,0,comm3d,ierr)
        if (iperiod >= nperiod) goto 60
    endif

    IF(mod(nn,ndump) == 0 .and. ipack >= 0) then 
        CALL DUMP(tcyc,outputfile)
    ENDIF

    if(tcyc > time_steady_pressure .and. tcyc-time_steady_pressure < 1.0d0*timestep &
        .and. ipack >= 0) then
        write(*,*)'STEADY LIMIT REACHED, DUMPING SOLUTION'
        CALL DUMP(tcyc,'stdysol')
    endif

    if(tcyc > time_oseen .and. tcyc-time_oseen < 1.01d0*timestep &
        .and. ipack >= 0) then
        write(*,*)'OSEEN LIMIT REACHED, DUMPING SOLUTION'
        CALL DUMP(tcyc,'oseesol')
    endif


    if (done_time) goto 60

end do  ! do nn=1,npts
!
!*******************************************************************************************
!                                 END OF THE COMPUTATION
!*******************************************************************************************

60 continue

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  total_timestep = MPI_WTIME() - total_timestep 

  ! PRF: Get the maximum time loop duration across all procs
  call MPI_ALLREDUCE(total_timestep,maxtimestep,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)

  ! PRF: The calculation time is execution time less I/O
  calctime = MPI_WTIME() - starttime - progtime

 

  if (npts.lt.0) goto 888
  
  temptime = MPI_WTIME()
  
888 continue

  if(myid.eq.0) then
     write(6,*)'alpha,beta = ',alpha,beta
  endif

  ! PRF: Finish time for program
  endtime = MPI_WTIME()

70 continue

  ! Close files
     close(9)
     close(29)
!     close(501)


!--modi0
  CALL PRINT_OUTPUT(tout) 



  CALL MPI_FINALIZE(ierr)
  
END PROGRAM ROCFIRE



!    ****************************************************************
!    Progress Report
!
!    This subroutine writes various error and maximum values to the 
!    screen during execution of the program
!
!    *****************************************************************  
SUBROUTINE PROGRESS(t,tout)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!---------------------------------------------------------------------
! Local variables
  INTEGER :: i,j,k,n,eqn 
  INTEGER :: iskip,iframe,writeflag,xs,zs 
  REAL*8 ::  t, tout, err1, diff, rmax, tgmax, tsmax,r1, xx, yy, zz
  REAL*8 ::  ap_area,total_area,amass_ap,amass_b
  REAL*8 ::  const,amass_total,equiv,heatflux
  REAL*8 ::  all_err1,allrmax,alltgmax,alltsmax,all_amass_b
  REAL*8 ::  tsmin,alltsmin,phinot,allphinot,phimax
  REAL*8 ::  all_ap_area,all_total_area,all_amass_ap, avg_speed, temptime
  REAL*8 ::  hx, hy, Tx, Ty, pressure, poo
  REAL*8 ::  hflux_ap, hflux_b,all_hflux_ap, all_hflux_b,hflux_total
  REAL*8 ::  vavg, Tavg, all_vavg, all_Tavg, mgas
  REAL*8 ::  heat_out,all_heat_out,heat_in,all_heat_in
  REAL*8 ::  rho_in,rho_out,all_rho_in,all_rho_out
  REAL*8 ::  temp_in,temp_out,all_temp_in,all_temp_out
  REAL*8 ::  kinetic_heat,beta_surf,c_h_g,adiabatic_temp
  REAL*8 ::  max_ft,min_ft,all_max_ft,all_min_ft
  INTEGER :: imax, jmax ,kmax, imax2, jmax2
!---------------------------------------------------------------------


  if(mod(ncyc,writemod) /=0 ) RETURN

  xs = drange(1)
  zs = drange(3)

! Calculate Maximum Values
  rmax = 0.0
  tsmin=200000.
  tgmax = 0.0
  tsmax = 0.0d0
  do j = drange(5), drange(6)
     do k = drange(3), drange(4)
        do i = drange(1), drange(2)
           r1   = rate(i,k,j,1)
           rmax = dmax1(r1,rmax)
!              tgmax = dmax1(f(i,k,j,1),tgmax)
           tsmax = dmax1(f(i,k,0,neqmax),tsmax)
           tsmin = dmin1(f(i,k,0,neqmax),tsmin)
           if(abs(f(i,k,j,1)) > abs(tgmax) ) then
              tgmax = abs(f(i,k,j,1))
              imax = i
              kmax = k
              jmax =j
           endif
        end do
     end do
  end do



  rmax = rmax*1.0d-6
! PRF: Get the max's and min's across all the processes, return them in all_XXXX
  if(nproc.gt.1) then
     call MPI_ALLREDUCE(rmax,allrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
     call MPI_ALLREDUCE(tgmax,alltgmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
     call MPI_ALLREDUCE(tsmax,alltsmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
     call MPI_ALLREDUCE(tsmin,alltsmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm3d,ierr)
     rmax = allrmax
     tsmax = alltsmax
  endif

  tgmax = alltgmax

!
! avg_speed is the period average of the regression speed, it is the data
! that it is usually compared against experiments, each processor will have a different
! value; we consider the one for processor (0)
!

  speed = phit(xs,zs)
  if(ipack == 1) then
     avg_speed = abs(phi(xs,zs)-phi_per)/(tout-t_per)
  else
     avg_speed = speed
  endif


! Calculate Total area, Total mass, and Equivelent
  eqn = neqmax
  hx = one/(2.0d0*dx)
  hy = - detady(0)/(2.0d0*dy)  !negative for solid
  hflux_ap = 0.0d0
  hflux_b = 0.0d0

  ap_area = 0.0d0
  total_area = 0.0d0
  amass_ap = 0.0d0
  amass_b = 0.0d0
  vavg = zero
  Tavg = zero
  heat_out = zero
  heat_in = zero
  rho_out = zero
  rho_in = zero
  temp_out = zero
  temp_in = zero
  max_ft = zero
  min_ft = 100.0d0
  do k=drange(3),drange(4)
     do i=drange(1),drange(2)
        const = sqrt(1.0+dphidx(i,k)**2+dphidz(i,k)**2)
        total_area = total_area+const
        Tx = (f(i+1,k,0,eqn)-f(i-1,k,0,eqn))*hx
        Ty = (-3.0d0*f(i,k,0,eqn)+4.0d0*f(i,k,1,eqn)-f(i,k,2,eqn))*hy
        vavg = vavg + q(i,k,nyv(3),3)
        Tavg = Tavg + f(i,k,nyv(0),1)
        heat_out = heat_out  + rhos(i,k,0) *phit(i,k)*f(i,k,0,neqmax)
        heat_in   =  heat_in + rhos(i,k,ny)*phit(i,k)*f(i,k,ny,neqmax)
        rho_out = rho_out + rhos(i,k,0)*phit(i,k)
        rho_in   =  rho_in + rhos(i,k,ny)*phit(i,k)
        temp_out = temp_out + f(i,k,0,neqmax)
        temp_in   =  temp_in + f(i,k,ny,neqmax)
        amass_ap = amass_ap+rhos(i,k,0)*phit(i,k)*xsurf(i,k,3)
        amass_b =  amass_b +rhos(i,k,0)*phit(i,k)*xsurf(i,k,2)
!---amass_ap = amass_ap+xsurf(i,k,3)
!---amass_b =  amass_b +xsurf(i,k,2)
        max_ft = max(max_ft,const*rb(i,k))
        min_ft = min(min_ft,const*rb(i,k))
        if (psi(i,k,0).ge.0.0d0) then
           ap_area=ap_area+const
           hflux_ap = hflux_ap + lambdas(i,k,0)/cp*(const**2*Ty - dphidx(i,k) * Tx)
        else
           hflux_b =  hflux_b  + lambdas(i,k,0)/cp*(const**2*Ty - dphidx(i,k) * Tx)
        endif
     enddo
  enddo

! PRF: Max,Mins across all processes
  if(nproc.gt.1) then
     call MPI_ALLREDUCE(total_area,all_total_area,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
     call MPI_ALLREDUCE(ap_area,all_ap_area,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
     call MPI_ALLREDUCE(amass_ap,all_amass_ap,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
     call MPI_ALLREDUCE(amass_b,all_amass_b,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
     call MPI_ALLREDUCE(hflux_ap,all_hflux_ap,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
     call MPI_ALLREDUCE(hflux_b,all_hflux_b,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
     call MPI_ALLREDUCE(vavg,all_vavg,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
     call MPI_ALLREDUCE(Tavg,all_Tavg,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
     call MPI_ALLREDUCE(heat_out,all_heat_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
     call MPI_ALLREDUCE(heat_in,all_heat_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
     call MPI_ALLREDUCE(rho_out,all_rho_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
     call MPI_ALLREDUCE(rho_in,all_rho_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
     call MPI_ALLREDUCE(temp_out,all_temp_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
     call MPI_ALLREDUCE(temp_in,all_temp_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
     call MPI_ALLREDUCE(max_ft,all_max_ft,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
     call MPI_ALLREDUCE(min_ft,all_min_ft,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm3d,ierr)
     total_area = all_total_area
     ap_area = all_ap_area
     amass_ap = all_amass_ap
     amass_b = all_amass_b
     hflux_ap = all_hflux_ap
     hflux_b = all_hflux_b
     vavg = all_vavg
     Tavg = all_Tavg
     heat_out = all_heat_out
     heat_in   = all_heat_in
     rho_out = all_rho_out
     rho_in   = all_rho_in
     temp_out = all_temp_out
     temp_in   = all_temp_in
     max_ft = all_max_ft
     min_ft = all_min_ft
  endif

! Some final calculations and adjustments
  ap_area=ap_area/(nx*nz)
  total_area=total_area/(nx*nz)
  amass_ap=amass_ap/(nx*nz)
  amass_b=amass_b/(nx*nz)
  amass_total=amass_ap+amass_b
  equiv=amass_b/amass_ap*beta
  hflux_ap = hflux_ap/(nx*nz)
  hflux_b = hflux_b/(nx*nz)
  hflux_total = hflux_ap+hflux_b
  vavg = vavg / (nx*nz)
  Tavg = Tavg / (nx*nz)
  heat_out = heat_out / (nx*nz)
  heat_in = heat_in / (nx*nz)
  rho_out = rho_out / (nx*nz)
  rho_in = rho_in / (nx*nz)
  temp_out = temp_out / (nx*nz)
  temp_in = temp_in / (nx*nz)
  mgas = pcyc*dim_fact/Tavg*vavg

  beta_surf = max(amass_ap/amass_b,beta)
  kinetic_heat = amass_ap/cp*(qg(1)+qheat_ap+qg(3)/beta_surf) &
       +  amass_b/cp*qheat_binder


  storage_term_2 = heat_out - heat_in - hflux_total

  c_h_g = - storage_term_2 + heat_in + kinetic_heat

  adiabatic_temp = (heat_in + kinetic_heat)/rho_out

  mass_flux = rho_out

  temptime = MPI_WTIME()


110 format(2x,19f14.6)
  poo = PRESSURE(tcyc)
  if(myid.eq.0) then
        write(6,601) TRIM(flag_comp),ncyc,tout,rb(0,0),Tavg,f(1,1,0,2),f(1,1,0,3),f(1,1,0,4),f(1,1,0,5),amass_total 

        if(ipack == 0) write(*,*)ncyc,'MAX-MIN SURF SPEED',&
             max_ft,min_ft,surfvel_max,surfvel_min
     
  endif
  if(myid == 0 .and. (mod(ncyc,frame) ==0) )then
        write(mflux,603) ncyc,tcyc,amass_total,ap_area,total_area,Tavg,vavg,mgas
  endif

601 format(a,2x,i6,18f12.6)
602 format(i4,1x,1p7e19.8,0p15f19.9)
603 format(i6,7(1x,',',0pf19.9))
604 format(3f14.6,8f14.6)
605 format(i8,4f12.6,a)


900 format(2x,i6,3f18.11)

  RETURN

END SUBROUTINE PROGRESS


!    ****************************************************************
!    Results Report
!    This subroutine writes solution results to files.
!    *****************************************************************    
SUBROUTINE PRINT_OUTPUT(tout)
  
  USE GLOBAL_DATA
  USE MYMPI
  
  IMPLICIT NONE

  ! Local variables
  INTEGER ::  i, j, k, iskip, eqn
  REAl*8 ::  tout, rmax, xx, yy, zz
  CHARACTER(LEN=21) :: fieldsfile,paramfile,zcfile,xcfile,ycfile,cfile,sfile
  
  ! PRF: For the time being, we are printing out only the X=0
  ! PRF: plane.  The rest of the procs must wait :P
  if(ddrange(1).eq.0) then
     WRITE(fieldsfile,'("fields",I3.3,".dat")') ddrange(3)
     WRITE(paramfile,'("param",I3.3,".dat")') ddrange(3)
     WRITE(zcfile,'("zcoord",I3.3,".dat")') ddrange(3)
     WRITE(xcfile,'("xcoord",I3.3,".dat")') ddrange(3)
     WRITE(ycfile,'("ycoord",I3.3,".dat")') ddrange(3)
     WRITE(cfile,'("coords",I3.3,".dat")') ddrange(3)
     WRITE(sfile,'("surface",I3.3,".dat")') ddrange(3)
     
900  format(2x,19f14.5)
901  format(2x,8f15.6)
902  format(2x,8f15.6)
     
     if (ipack.eq.1) iskip=2
     iskip=1
     
     ! Write out x & y physical coordinates
!     open (UNIT=93,FILE=cfile,STATUS='UNKNOWN')
!     do j=drange(5),drange(6)
!        do k=drange(3),drange(4),iskip
!           do i=drange(1),drange(2),iskip
!              xx = x(i)
!              zz = z(k)
!              yy = y(j) + phi(i,k)
!              write(93,900) xx, zz, yy
!           end do
!        end do
!     end do
!     close(93)
!
     ! Write out coordinates
!     open (UNIT=91,FILE=xcfile,STATUS='UNKNOWN')
!     do i = drange(1), drange(2)
!        write(91,901) i, x(i)
!     end do
!     close(91)
!
!     open (UNIT=90,FILE=zcfile,STATUS='UNKNOWN')
!     do k = drange(3), drange(4)
!        write(90,901) k, z(k)
!     end do
!     close(90)
!
!     open (UNIT=92,FILE=ycfile,STATUS='UNKNOWN')
!     do j = drange(5), drange(6)
!        write(92,901) j, y(j), detady(j), deta2dy2(j)
!     end do
!     close(92)
!     
     ! Write out the surface data...
!     open (UNIT=95,FILE=sfile,STATUS='UNKNOWN')
!     do k = drange(3), drange(4)
!        do i = drange(1), drange(2)
!           write(95,901) x(i), z(k), rb(i,k), phi(i,k),& 
!                dphidx(i,k),dphi2dx2(i,k),psi(i,k,0),&
!                f(i,k,0,neqmax),vvel(i,k,0)
!        end do
!     end do
!     close(95)

     ! Write out field variable data...
!     open (UNIT=10,FILE=fieldsfile,STATUS='UNKNOWN')
!     do j = drange(5), drange(6)
!        do k = drange(3), drange(4), iskip
!           do i = drange(1), drange(2), iskip
!              xx = x(i)
!              zz = z(k)
!              yy = y(j) + phi(i,k)
!              rmax = rate(i,k,j,1)*1.0d-6
!              write(10,900) xx, zz, yy, &
!                   (f(i,k,j,eqn), eqn=1,neqmax),&
!                   rmax, vvel(i,k,j),&
!                   phi(i,k), psi(i,k,j), rb(i,k), tout
!           end do
!        end do
!     end do
!     close(10)
!
!-modi0
     write(6,*) 'neqmax', neqmax

     ! Write out field variable data...
     open (UNIT=10,FILE=fieldsfile,STATUS='UNKNOWN')


      write (10,*)'Variables = "x" "y" "f1" "f2" "f3" "f4" "f5" "rmax" "vvel" "phi" "psi" "rb" "tout"'

      write(10,*) 'zone  I= ', drange(2)-drange(1)+1,'  J=', drange(6)-drange(5)+1
     k=0 
     do j = drange(5), drange(6)
           do i = drange(1), drange(2), iskip
              xx = x(i)
              yy = y(j) + phi(i,k)
              rmax = rate(i,k,j,1)*1.0d-6
              write(10,900) xx, yy, &
                   (f(i,k,j,eqn), eqn=1,neqmax),&
                   rmax, vvel(i,k,j),&
                   phi(i,k), psi(i,k,j), rb(i,k), tout   
           end do
        end do
!     end do
     close(10) 



  endif

  ! PRF: This just keeps the output looking nice, and each process 
  ! PRF: finishes at the same time.
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN

END SUBROUTINE PRINT_OUTPUT

!    ****************************************************************
!    Variable Input Routine
!    Written by Gregory M. Knott - 22 January 2000

!    This subroutine takes the inputs required for the modeling code
!    from the file "input.txt"
!    *****************************************************************

SUBROUTINE INPUTS

USE GLOBAL_DATA
USE MYMPI

IMPLICIT NONE

! Local variables
integer :: inunit,i,j
CHARACTER :: comt*100
REAL*8 ::  ru,lewis_y,mw


inunit = 120
open(inunit,FILE='input.txt',STATUS='old')

! Read main inputs
! General Code Parameters
read(inunit,*) comt
read(inunit,*) comt
read(inunit,*) comt
read(inunit,*) nyprocs
read(inunit,*) irestart
read(inunit,*) restartfile
read(inunit,*) outputfile
read(inunit,*) npts
read(inunit,*) timestep
read(inunit,*) tstop
read(inunit,*) dflag
read(inunit,*) writemod
read(inunit,*) ndump
read(inunit,*) frame
read(inunit,*) n_print_terms
read(inunit,*) ipack
read(inunit,*) ncyc_oseen
read(inunit,*) nperiod
read(inunit,*) period_begin
read(inunit,*) ncyc_steady
read(inunit,*) ncyc_constant
!
!    Grid Parameters
!    note the dimensions in x and z are the same
!
read(inunit,*) comt
read(inunit,*) c1y
read(inunit,*) nnxx
read(inunit,*) nnyy
read(inunit,*) lvlx
read(inunit,*) lvly
read(inunit,*) xend
read(inunit,*) yend

! Pressure and Pressure Disturbance Parameters
read(inunit,*) comt
read(inunit,*) press
read(inunit,*) epi_p
read(inunit,*) omg_p
read(inunit,*) phase_p

if(mod(phase_p,one) > 1.d-8) then
    write(*,*)'ERROR OLD INPUT DECK'
    STOP
endif

write(*,*)'NEW VARIABLE, PHASE_P',phase_p

! Amount AP and Homginization parameter
read(inunit,*) comt
read(inunit,*) AP_frac
read(inunit,*) alpha
read(inunit,*) bind_sizes
read(inunit,*) (bind_diam(j),j=1,bind_sizes)
read(inunit,*) (bind_mass(j),j=1,bind_sizes)
read(inunit,*) premixed

! Solid Phase Physical Parameters
read(inunit,*) comt
read(inunit,*) tcold
read(inunit,*) rho_ap
read(inunit,*) rho_binder

! Time delays
read(inunit,*) comt
read(inunit,*) time_steady_pressure
read(inunit,*) time_oseen

! Global Gas Phase Properties
read(inunit,*) comt     
read(inunit,*) cp       !Gas and Solid heat capacity
read(inunit,*) pr       !Prandtl Number
read(inunit,*) ru       !Gas Constant
read(inunit,*) mw       !Average gas MW

! Close file
close(inunit)

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!............................................................................
! CALCULATE RELEVANT PARAMETERS
!............................................................................

dim_fact = ru/mw * 4186.8d0/101325.0d0*1000.0d0 ! kcal/m^3->atm
dim_fact=1.0d0/dim_fact  !Dimensional factor

pi = ACOS(-1.0d0)

! Convert microns to cm
xend = xend/1.0d4
yend = yend/1.0d4
zend = xend

!    Set-up the Lewis numbers
lewis=1.0d0;lewis(2)=lewis_y
!
!this is the accelration of gravity which should be set to zero if not includede
!
g_accel = zero


!----------------------------------------------------
RETURN
END SUBROUTINE INPUTS
!********************************************************
