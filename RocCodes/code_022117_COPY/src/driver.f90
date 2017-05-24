! ********************************************************************
!                          ROCFIRE3D START POINT 
! ********************************************************************
! ********************************************************************
!
! Heterogeneous Propellant 3D Combustion Code
!
! This is a legacy code originally developed at the
! DOE-ASCI Center for Simulation of Advanced Rockets
! University of Illinois at Urbana-Champaign
!
! Author: Thomas L. Jackson
! Date last Modified: 16 January 2001
! Parallelized by M. Campbell on May 9, 2001
! Extensively worked on by Luca Massa
!
! ////////////////////////////////////////////////////////////////////


! comments that may help when reading code:
!     f(i,k,j,1)  = T    (gas phase temperature)
!     f(i,k,j,2)  = Y_O  (gas phase oxidizer X)
!     f(i,k,j,3)  = Y_Z  (gas phase intermediate Z)
!     f(i,k,j,4)  = Y_F  (gas phase fuel Y)
!     f(i,k,j,5)  = T    (solid phase temperature)
!
!     MQ*** = ???
!
!     turn off radiation by making NCYCRAD = 1 gazillion ???
!
!     G_*** is variables for printing, I think (Gathered from all threads)
!
!     poo has something to do with pressure


PROGRAM ROCFIRE  
  USE GLOBAL_DATA
  USE BC_DATA
  USE MYMPI
  use statistics
  Use restart_dump
  Use Timing
  use implicit
  USE radiation

  IMPLICIT NONE

  include "parallel.h"

!--------------------------------------------------------------------
! Local variables
  INTEGER ::  filenum, i, i5, m, n,j,k,idump, nn
  INTEGER ::  ii,imax,eqn,ncyc_steady_ipack
  INTEGER ::  sum1,sum2,sum_K,npts_run,myploc
  REAL*8  ::  tout, t,yy
  REAL*8  ::  error, u3, v3, w3
  REAL*8  ::  per_test
  REAL*8  ::  vol,vol_tot,vol_ap,vol_bind,vol_pack
  REAL*8  ::  ratio
  REAL*8  ::  ddlx,ddly,ddlz
  REAL*8  ::  rvtot_K,rvexp_K,beta_K,diff_beta_K
  REAL*8  ::  th,lb1,lb3,coe_rho
  REAL*8  ::  ymax,zmax,xmax,radmax
!
! PRF: Timing variables
!
  REAL*8  ::  endtime,inputtime,ictime,quad2time
  REAL*8  ::  packtime,setuptime,gridtime,udpsitime,iudpsitime
  REAL*8  ::  spherebctime,quadtime,calctime,optime
  REAL*8  ::  txytime, phitime,bctime,restarttime
  REAL*8  ::  progtime,temptime,psitime,ictime1,maxsetuptime
  REAL*8  ::  outtime,decomptime,total_timestep,maxtimestep,pressure

  LOGICAL :: MaxSimuTimeReached,flag1,flag2

! PRF: Process specific output files
  CHARACTER(LEN=21) :: rstartfile,statfile,headerfile,timingfile
  CHARACTER(LEN=21) :: moviefile
!--------------------------------------------------------------------------

  call MPI_INIT(ierr)

  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Prid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  starttime = MPI_WTIME()
  done_time = starttime
  maxtime = 1d99
  steadyconvergence = .false.
  numerdivergence = .false.

  dooutput = 1

  Nimposed = 2
!
! Read in from filenames.list the example file to run
!
  inreadfilename = 601
  open(600,file='filenames.list',ACTION='READ')
  read(600,*) inputfile
  do while (scan(inputfile,'!')==1)
     read(600,*) inputfile
  enddo
  if (myid==0) write(6,*)'In driver, inputfile=', TRIM(inputfile)
  open(inreadfilename,file=TRIM(inputfile),ACTION='READ')

  read(inreadfilename,*) inputfile
  inputfile_name=TRIM(inputfile)
  if (myid==0) write(6,*)'In INPUTS, inputfile=', inputfile_name

  read(inreadfilename,*) inputfile
  chem_inputfile_name=TRIM(inputfile)
  if (myid==0) write(6,*)'In CHEMICAL, inputfile=', chem_inputfile_name

  read(inreadfilename,*) inputfile
  pack_inputfile_name=TRIM(inputfile)
  if (myid==0) write(6,*)'In PACK, inputfile=', pack_inputfile_name
!
!
  CALL INPUTS
!
! Set maximum size
!
  neqmax = maxsize
  neqgas = neqmax-1
  nspecies = neqgas-1
  irxn = nspecies-1
  if(NORADIATION) irxn = nspecies

  if(myid == 0) then
     write(6,*)'maxsize=',maxsize
     write(6,*)'neqmax=',neqmax
     write(6,*)'neqgas=',neqgas
     write(6,*)'nspecies=',nspecies
     write(6,*)'irxn=',irxn
  endif

!skip flags, should be all false
  is2D = .FALSE.
  is_firstorder =        .true.
  skipTVD =                    .true.
  skipDeformation =      .not. .true.
  skipDeformationFluid = .not. .true.
  skipTemperature =      .not. .true.
  skipViscDeriv =        .not. .true.
  skipContSource =       .not. .true.  .or. skipTemperature
  skipPressure =         .not. .true.
  skipVarDensity =       .not. .true.
  skipCrossflow =        .true.
  skipUnsteady  =        .not. .true.
  skipNorthNeumanBcTemp =.not. .true. .and. .not. skipCrossflow

! reset some flags, TLJ
  is_firstorder = .true.
  skipTVD = .true.
  if (NORADIATION) skipTVD = .false.

  InitialPressure = press
  initialTimestep = timestep


  call read_chem_info

  IF (ipack <= 0) THEN
     dlength = xend
     dly = dlength
     period = dlength
     if(iwflag /= 1) then
        alpha =  alpha/(rho_binder/rho_ap*(1.0d0-alpha)+alpha)
        iwflag = 1
     end if
     ncyc_steady_ipack = ncyc_steady
     if(ipack == 0) then 
        coe_d = zero
        is2D = .true.
        issymmetric = .TRUE.
        write(6,*)'Calling SAND_PACK'
        CALL SAND_PACK (packtime)
     elseif (ipack < 0) then 
        coe_d = merge(zero,one,ipack == -1)      !steady solid vs unsteady solid
        ipack = -1
        ncyc_oseen = npts+1
        nnxx = 1
        lvlx = 1
        writemod = 250
        irestart = 0
        ndump = npts+1
        if(nproc >1)then
           write(*,*)'CASE ipack = -1 runs in single proc mode only'
           stop 'driver.f90'
        endif
     endif
  ELSE
     coe_d = one     
     issymmetric = .false.
     if (ndim==2) then
        is2D = .true.
        CALL READ_DISCS (packtime)
     endif
     if (ndim==3) then
        is2D = .false.
        CALL READ_3DPACK (packtime)
     endif
     if(.not. skipCrossFlow) call READ_MATRIX
  ENDIF

  call global_parameter
  call RADIATION_hardwire

! this call finds the lambda of the blend
  CALL BLEND((/ts_AP_H,ts_BD_H,ts_alum_H/),(/lambda_ap,lambda_binder,lambda_al/),&
       &(/rho_ap,rho_binder,rho_al/),lambda_eff,rho_eff,3)


  if (myid==0) then
     write(6,*)'homog=',ts_AP_H,ts_BD_H,ts_alum_H
     write(6,*)'lambda=',lambda_ap,lambda_binder,lambda_al
     write(6,*)'blend=',lambda_eff,rho_eff
     write(6,*)'alp_V=',alp_V(1:3)
     write(6,*)'alp_W=',alp_W(1:3)
  endif


!
! ASSIGN number of processor in each dimension: pdims(1:3)    
!
  CALL GET_PROC_DECOMP(1)
  CALL GET_MG_DECOMP
  CALL GET_PROC_DECOMP(2)

  CALL ALLOCATE_MAIN

! Assign ouput file names
  flagfile=press

  WRITE(statfile,'("stats",I3.3,".dat")') flagfile
  WRITE(errfile,'("history",I3.3,".dat")') flagfile
  WRITE(timingfile,'("timing",I3.3,".dat")') flagfile


!
! Setup Grid and Initial Guess
!
  CALL GRIDSETUP

  if(ipack == 0) then
     iloc=  0
     myploc = 0
     DO i = drange(1), drange(2)
        if(x(i) > xloc .AND. x(i) - xloc < dx ) then
           iloc = i
           myploc = myid
        endif
     END DO
     CALL MPI_ALLREDUCE(myploc,ploc,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
  endif
  10 format("% cycle   time (s)    rb (cm/s)  rb_avg_200   rb_avg_1k   &
      &rb_avg_5k  rb_avg_25k  rb_avg_100k  rb_avg_1M   Tsurf_avg   &
      &phi_avg     Tsurf_max   Tsurf_min   Tgas_max   &
      &equiv      rmax      surf_area     mass_flux")
  IF (irestart == 0) THEN
     if (myid == 0) write(6,*)'Calling INIT_FROM_SCRATCH'
     CALL INIT_FROM_SCRATCH
     if (myid == 0) write(6,*)'End INIT_FROM_SCRATCH'
     if(myid == 0) then
       open(UNIT=9,FILE=errfile,STATUS='UNKNOWN',POSITION='REWIND')
       if( ipack==1 ) write(9,fmt=10)
       endif
  ELSEIF (irestart == 1) THEN
     CALL RESTART(restartfile)
     if(myid == 0) open(UNIT=9,FILE=errfile,STATUS='UNKNOWN',POSITION='APPEND')
     if(myid == 0) write(9,*) '%restarted ncyc= ',ncyc,'Time = ',tinit
  ELSEIF (irestart == -1) THEN
     CALL RESTART_1D
     if(myid == 0) then
       open(UNIT=9,FILE=errfile,STATUS='UNKNOWN',POSITION='REWIND')
       if( ipack==1 ) write(9,fmt=10)
       endif
  ENDIF


  CALL ALlOCATE_SOLVER   !the allocation were split up to allow for visualization
  CALL GET_TIME_COEFF
  CALL GRID_SPEED

  if(Nprobes > 0.and. ipack == 1) then
     CALL PROBESETUP
     if(irestart == 1) then
        call PROBERESTART
     endif
  endif

  CALL PRN_GLOB_VARS

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  total_timestep = MPI_WTIME()

! PRF: Get the maximum setuptime across all procs
  CALL MPI_ALLREDUCE(setuptime,maxsetuptime,1,&
     &MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)

!*********************************************************************
!*********************************************************************
!                              START COMPUTATION

  gausstime = zero
  commtime = zero
  startVIZ2 = .TRUE.
  t = tinit
  tout = t
  tcyc = t   
  CALL CONTINUITYSOURCEINIT
  !CALL FourierExtVelUpdate(0,1)


  npts_run = npts - ncyc + 1
  ncyc_P = 0
  NCYCRAD = 300
  ncyc_init = 500
  pi = acos(-1.0d0)
  !if(myid==0) print *, '\n\n----- begin time step ------'


  TimeIntegrationLoop : do nn=1,npts_run

     call Taketime(1000000)

     ! shear only works for Oseen at the moment
     shear_rate = shearRateAmp
     do j = drange(5), drange(6)
        shear(j) = shearRateAmp*y(j)*&
               & SIN(ncyc*timestep*shearFreq*2.*pi)
     enddo

     if(ipack == -1 )then
        ncyc_P = ncyc_P+1
        if(mod(ncyc_P,2600) == 0 .and. ncyc > 1000) then   !cancel
           ncyc_P = 0
           if(PRESSVARY) then
              if(press < 7d1) then
                 press = press + 5d0   
              else
                 press = InitialPressure
                 KrierFACTOR = KrierFACTOR/2d0
              end if
           elseif(HEATVARY) then
              KrierFACTOR = KrierFACTOR/2d0
           end if
        end if
     elseif(ipack == 0 .and. PRESSVARY .or. HEATVARY )then
        ncyc_P = ncyc_P+1
        if(mod(ncyc_P,365) == 0 .and. ncyc > 900) then   !cancel
           ncyc_P = 0
!!>           timestep = initialTimestep*2d0
           if(PRESSVARY) then
              press = press + 1d1
           elseif(HEATVARY) then
              KrierFACTOR = KrierFACTOR/2d0
           end if
        end if
     end if

     if(mod(nn-1,1) == 0) saved_LU_vars = .FALSE.

     dt = timestep

     if(ipack <= 0) coe_d = merge(zero,one, nn > ncyc_steady_ipack)

     ncyc=ncyc+1
     CALL GET_TIME_COEFF

     ncyc_run = nn
     pcyc = pressure(tcyc(2))  !spatially constant


     CALL OLDSOLUTION
     CALL BCQ(1)
     call update_phi(zero,zero,0) !update only the derivatives

!    Update radiation field if included
     if (NORADIATION) BCradheat = 0d0
     if (NORADIATION.eqv..FALSE.) then
     if(ncyc >= NCYCRAD) then
        BCradheat = 0d0
        CALL CURB_TXY
        call LAMBDA(0)
        CALL RADIATION_SOURCE
     else
        BCradheat = 5d1
     end if
     end if

!--- Update equations
     CALL BC
     CALL UPDATE_TXY  ! calls gmres and TXY_UPDATE (seems to apply the derivative to update f)

     CALL CURB_TXY    ! seems to put limits on f

!    Update Boundary Conditions
     CALL BC
     CALL PYROLYSIS
     CALL BCQ(1)
     call VELOCITY(0)   ! compute conserved velocities


     CALL CONTINUITYSOURCE

!--- Update Mid point vars
     CALL MidPointDensityConductivity


!--- Update surface
     tcyc(1) = t
     CALL PROPAGATE_SURF(t,tout,merge(0,merge(-1,1,skipTemperature),&
          skipUnsteady),MaxSimuTimeReached)
     tcyc(2) = tout
     WRITE(timeString,'(1PE11.5)') tcyc(1)

     !CALL FourierExtVelUpdate(0,1)
     !CALL FourierExtVelUpdate(1,1)
     if(doprojection) then
        CALL BCQ(1)
        CALL UPDATE_QXY
        CALL PROJECTION(0)
     endif

     CALL BCQ(2)
     !CALL FourierExtVelUpdate(5,1)

     !!$     call Test_MGcoeffs
     if (mod(ncyc,2*writemod) == 0) CALL PROJECTION(1)

     if(ipack == 1) CALL Stat_Average
     CALL PROGRESS(t,tout)


     if (ipack.eq.1) then
        iperiod = int(phi(0,0)/period) - int(phi_begin/period)
        call MPI_BCAST(iperiod,1,MPI_INTEGER,0,comm3d,ierr)
        if (iperiod >= nperiod) goto 60
     endif

     !if(ipack == 0) then
     !   CALL vizSAND(mod(ncyc_run-1,frame)==0)
     !else
     !   CALL vizPACK3D(mod(ncyc_run-1,frame)==0)
     !end if

     if(ipack == 1 .and. ( mod(ncyc_run-1,frame)==0 .or. ncyc == npts ) ) then
       ! call shaneWriteField( ncyc_run )
       call shaneWriteSlice( ncyc_run )
     end if

     !IF(mod(ncyc_run,ndump) == 0) then 
     !   CALL DUMP(tout,outputfile)
     !ENDIF

     call TimeReport

     call alldone

     call MPI_ALLREDUCE(steadyconvergence,flag1,1,MPI_LOGICAL,MPI_LOR,comm3d,ierr)
     call MPI_ALLREDUCE(numerdivergence,  flag2,1,MPI_LOGICAL,MPI_LOR,comm3d,ierr)

     !if(myid==0) print *, '-----  end time step  ------'

     if (MaxSimuTimeReached.or.flag1.or.flag2) then
        if(myid == 0) print*,'EXITING','MaxSimuTimeReached',MaxSimuTimeReached,&
             &'flag1',flag1,'flag2',flag2
        if(ipack == -1) then
        else
           goto 60
        end if
     endif

     if(ncyc == npts)  goto 60

  end do TimeIntegrationLoop  ! do nn=1,npts
!
!*********************************************************************
!                                 END OF THE COMPUTATION
!*********************************************************************

60 continue

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  total_timestep = MPI_WTIME() - total_timestep 

! PRF: Get the maximum time loop duration across all procs
  call MPI_ALLREDUCE(total_timestep,maxtimestep,1,&
      &MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)

! PRF: The calculation time is execution time less I/O
  calctime = MPI_WTIME() - starttime - progtime

! PRF: Output notification
  if(dooutput.eq.1) then
     if(myid.eq.0) then
        write(*,*)
        write(*,*) 'Writing output files...'
     endif

! Print out values at final time step
     temptime = MPI_WTIME()
     CALL PRINT_OUTPUT(tout)
     temptime = MPI_WTIME() - temptime
     outtime = temptime
  endif

  if (npts.lt.0) goto 888

  temptime = MPI_WTIME()

888 continue

  if(myid.eq.0) then
     write(6,*)'alpha,beta,T_melt = ',alpha,beta,T_melt
  endif

! PRF: Finish time for program
  endtime = MPI_WTIME()

  if (myid.eq.0) then
     open (UNIT=87,FILE=timingfile,STATUS='UNKNOWN')       
     write(87,*)
     write(87,*) "--------------------------------------------------------------"
     write(87,*) "INITIALIZATION TIMES:"
     write(87,*) "   Input                :=",inputtime ,'seconds'
     write(87,*) "   PackTime             :=",packtime,'seconds'
     write(87,*) "   IC-I                 :=",ictime1,'seconds'
     write(87,*) "   IC-II                :=",ictime,'seconds'
     write(87,*) "   Quadrant             :=",quadtime,'seconds'
     write(87,*) "   DecompTime           :=",decomptime,'seconds'
     write(87,*) "   GridTime             :=",gridtime,'seconds'
     write(87,*) "   InitPsiTime          :=",iudpsitime,'seconds'
     if(doparmdump.eq.1) then
        write(87,*) "   OutputTime           :=",optime,'seconds'
     endif
     write(87,*) "   TotalSetupTime       :=",setuptime,'seconds'
     write(87,*) "   MaxSetupTime         :=",maxsetuptime,'seconds'
     write(87,*) "-------------------------------------------------------------"
     write(87,*) "LOOP TIMES:"
     write(87,*) "   TXY                  :=",txytime,'seconds'
     write(87,*) "   Psi                  :=",psitime,'seconds'
     write(87,*) "   Phi                  :=",phitime,'seconds'
     write(87,*) "   BCond                :=",bctime,'seconds'
     if(doprogress.eq.1) then
        write(87,*) "   Progress             :=",progtime,'seconds'
     endif
     write(87,*) "   Communication time   :=",commtime(1),'seconds'
     write(87,*) "   Gausscomm  time      :=",gausstime(1),'seconds'
     write(87,*) "   TotalTimeLoop        :=",total_timestep,'seconds'
     write(87,*) "   MaxTotalTimeLoop     :=",maxtimestep,'seconds'
     write(87,*) "------------------------------------------------------------"
     if(dooutput.eq.1) then
        write(87,*) "OUTPUT TIME:"
        write(87,*) "   FinalOutput               :=",outtime,'seconds'
        write(87,*) "---------------------------------------------------------"
     endif
     write(87,*) "   Total Execution Time :=",(endtime-starttime),'seconds'
     write(87,*) "------------------------------------------------------------"     
     close (87)
  endif

70 continue

! Close files
  if(myid.eq.0.and.doprogress.eq.1) then
     close(9)
!     close(501)
  endif

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
  use statistics

  IMPLICIT NONE

!---------------------------------------------------------------------
! Local variables
  INTEGER :: i,j,k,n,eqn ,l
  INTEGER :: iskip,iframe,writeflag,xs,zs ,xe,ze, ys,ye
  REAL*8 ::  t, tout, err1, diff, rmax, tgmax, tsmax,r1, xx, yy, zz
  REAL*8 ::  ap_area,total_area,amass_ap,amass_b
  REAL*8 ::  const,amass_total,equiv,heatflux
  REAL*8 ::  all_err1,allrmax,alltgmax,alltsmax,all_amass_b
  REAL*8 ::  tsmin,alltsmin,phinot,allphinot,phimax
  REAL*8 ::  all_ap_area,all_total_area,all_amass_ap, temptime
  REAL*8 ::  hx, hy, Tx, Ty, pressure, poo
  REAL*8 ::  hflux_ap, hflux_b,all_hflux_ap, all_hflux_b,hflux_total
  REAL*8 ::  vavg, Tavg, all_vavg, all_Tavg, mgas
  REAL*8 ::  heat_out,all_heat_out,heat_in,all_heat_in
  REAL*8 ::  rho_in,rho_out,all_rho_in,all_rho_out
  REAL*8 ::  temp_in,temp_out,all_temp_in,all_temp_out
  REAL*8 ::  kinetic_heat,beta_surf,c_h_g,adiabatic_temp
  REAL*8 ::  max_ft,min_ft,all_max_ft,all_min_ft
  REAL*8 ::  tommax,tomavg,tomterm,prandtl_cp,vcomm,Ux,Uy
  REAL*8 ::  MaxTemp,MinTemp,maxTempgas
  INTEGER :: imax, jmax ,kmax, imax2, jmax2,imin, jmin ,kmin
  character*30 :: fname
!---------------------------------------------------------------------
  call MPI_ALLREDUCE(rmax,allrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)


  if(mod(ncyc_run-1,writemod) /=0 ) RETURN

  prandtl_cp = pr/cp

  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  ys = drange(5)
  ye = drange(6)

! Calculate Maximum Values:
!   total heat release, rmax
!   min surface temp, tsmin
!   max surface temp, tsmax
!   max gas-phase temp, tgmax
  rmax  = 0.0
  tsmin = 200000.0
  tsmax = 0.0d0
  tgmax = 0.0d0
  do k = drange(3), drange(4)
     do i = drange(1), drange(2)
        j = 0
        tsmax = dmax1(f(i,k,0,neqmax),tsmax)
        tsmin = dmin1(f(i,k,0,neqmax),tsmin)
        do j = drange(5), drange(6)
           r1   = sum(yf(1,1:irxn)*rate_out(i,k,j,1:irxn))
           rmax = dmax1(r1,rmax)
           tgmax = max(f(i,k,j,1),tgmax)
        end do
     end do
  end do
  rmax = rmax*1.0d-6

! PRF: Get the max's and min's across all the processes, return them in all_XXXX
  if(nproc.gt.1) then
     call MPI_ALLREDUCE(tgmax,alltgmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
     call MPI_ALLREDUCE(tsmax,alltsmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
     call MPI_ALLREDUCE(tsmin,alltsmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm3d,ierr)
     call MPI_ALLREDUCE(rmax,allrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
     tgmax = alltgmax
     tsmax = alltsmax
     tsmin = alltsmin
     rmax = allrmax
  endif
!
! speed is the regression speed, it is the data
! that it is usually compared against experiments, 
! each processor will have a different
! value; we consider the one for processor (0)
!
  speed = phit(xs,zs)

! Calculate Total area, Total mass, and Equivelent
  eqn = neqmax
  hx = one/(2.0d0*dx)
  hy = - detady(0)/(2.0d0*dy)  !negative for solid
  hflux_ap = 0.0d0
  hflux_b = 0.0d0

  total_area = 0.0d0
  ap_area = 0.0d0
  amass_ap = 0.0d0
  amass_b = 0.0d0
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

     ! compute total integrated surface area
     total_area = total_area+const

     heat_out = heat_out  + rhos(i,k,0) *phit(i,k)*f(i,k,0,neqmax)
     heat_in   =  heat_in + rhos(i,k,ny)*phit(i,k)*f(i,k,ny,neqmax)
     rho_out = rho_out + rhos(i,k,0)*phit(i,k)
     rho_in   =  rho_in + rhos(i,k,ny)*phit(i,k)
     temp_out = temp_out + f(i,k,0,neqmax)
     temp_in   =  temp_in + f(i,k,ny,neqmax)

     ! equivalence ratio = surface integrated flux of fuel divided
     ! by surface integrated flux of AP, normalized by the
     ! stoichiometric proportion of these fluxes. Can be called
     ! the integrated flux based equivalence ratio.
     ! Note: this needs to be updated if xsurf in solid_prop.f90
     ! is changed.
     amass_ap = amass_ap+rhos(i,k,0)*phit(i,k)*xsurf(i,k,2)
     amass_b =  amass_b +rhos(i,k,0)*phit(i,k)*xsurf(i,k,4)

     ! Compute max and min of df/dt
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
  total_area=total_area/(nx*nz)
  ap_area=ap_area/(nx*nz)
  amass_ap=amass_ap/(nx*nz)
  amass_b=amass_b/(nx*nz)
  amass_total=(amass_ap+amass_b)/stat_rb(1)
  equiv=amass_b/amass_ap*beta
  hflux_ap = hflux_ap/(nx*nz)
  hflux_b = hflux_b/(nx*nz)
  hflux_total = hflux_ap+hflux_b
  heat_out = heat_out / (nx*nz)
  heat_in = heat_in / (nx*nz)
  rho_out = rho_out / (nx*nz)
  rho_in = rho_in / (nx*nz)
  temp_out = temp_out / (nx*nz)
  temp_in = temp_in / (nx*nz)


  temptime = MPI_WTIME()


110 format(2x,19f14.6)
  poo = PRESSURE(tcyc(2))   
  call MPI_ALLREDUCE(maxval(f(xs:xe,zs:ze,0,1)),maxTemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
  call MPI_ALLREDUCE(minval(f(xs:xe,zs:ze,0,1)),minTemp,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm3d,ierr)

  if(myid.eq.0) then

     if(abs(ipack) == 1) then

        if(ipack == 1) then
           write(6,601) TRIM(flag_comp),ncyc,tout,stat_rb(1:2)/dble(ncyc_stat(1:2)),stat_Ts(1) &
                ,pcyc,phiAvg,maxTemp,minTemp,tgmax,equiv,rmax,total_area,amass_total
           write(9,'(i7,1p123e12.4)') ncyc,tout,stat_rb(1:nstat)/dble(ncyc_stat(1:nstat)),&
                &stat_Ts(1)/dble(ncyc_stat(1)),phiAvg,maxTemp,minTemp,tgmax,equiv,&
                &rmax,total_area,amass_total
        else
           write(6,'(i8,1p123e14.4)')ncyc,tout,rb(0,0),f(0,0,0,1),f(0,0,ny,1),press,alp_W,coe_d
        end if

        if(mod(ncyc,4*writemod) <=1 .and. ipack==-1) then
           write(fname,'(a,"_",1pe7.1,"_",1pe7.1,a)')trim(scase),KrierFACTOR,press,'.dat'
           close(1233)
           open(1233,file = trim(fname))
           write(1233,'(/"% ",4(a12))')'Press','KrierFactor','Burnrate','Tflame'
           write(1233,'(2x,4(1pe12.4))')press,KrierFACTOR,rb(0,0),f(0,0,ny,1)
           close(1233)
        end if

     else

        if(ncyc ==1 .or. mod(ncyc,10*writemod) == 0) then
           write(6,'(a8,44a11)') 'Ncyc','Time','Maxres(1)','Maxres(2)',&
             &'Maxres(3)','T_g.max','T_s,max','phi','Speed','rb','Max speed',&
             &'Min Speed'
        end if
        write(6,701) ncyc, tout, maxres(1:3), tgmax, tsmax, &
             phi(xs,zs), speed, rb(xs,zs), surfvel_max, surfvel_min
        write(9,702) tout,temptime - starttime,log10(max(maxres(1),maxres(2))),&
             tgmax,tsmax,phi(xs,zs), speed, rb(xs,zs), surfvel_max, surfvel_min
        if(mod(ncyc,50*writemod) == 0) then
           write(6,'(/,i5,4a12,/,5x,1p4e12.4,/)')ncyc,'XLOC','COE_D',&
      &'COE_DT','Press',xloc,coe_d,coe_dt,Press
        end if
     endif

  endif


601 format(a,2x,i6,30f14.6)
602 format(i4,1x,1p7e19.8,0p15f19.9)
604 format(16HTom value is -- ,3f19.9)
605 format(i8,4f12.6,a)
701 format(i8,1p44e11.3)
702 format(1p3e12.4,0p15f15.7)

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
  INTEGER ::  i, j, k, iskip, eqn, m
  REAL*8 ::  tout, rmax, xx, yy, zz
  CHARACTER(LEN=21) :: fieldsfile,paramfile,zcfile,xcfile,ycfile,cfile,sfile
  REAL*8  :: poo
  REAL*8  :: pre(maxsize)
  CHARACTER(LEN=50) :: ratefile


! PRF: For the time being, we are printing out only the X=0
! PRF: plane.  The rest of the procs must wait :P
  if(dSerialRange(1).eq.0) then
     WRITE(fieldsfile,'("fields",I3.3,".dat")') dSerialRange(3)
     WRITE(paramfile,'("param",I3.3,".dat")') dSerialRange(3)
     WRITE(zcfile,'("zcoord",I3.3,".dat")') dSerialRange(3)
     WRITE(xcfile,'("xcoord",I3.3,".dat")') dSerialRange(3)
     WRITE(ycfile,'("ycoord",I3.3,".dat")') dSerialRange(3)
     WRITE(cfile,'("coords",I3.3,".dat")') dSerialRange(3)
     WRITE(sfile,'("surface",I3.3,".dat")') dSerialRange(3)

900  format(2x,19f14.5)
901  format(2x,8f15.6)
902  format(2x,8f15.6)

     if (ipack.eq.1) iskip=2
     iskip=1

! Write out physical and grid parameters
     open (UNIT=89,FILE=paramfile,STATUS='UNKNOWN')     
     write(89,*) tout
     write(89,*) nx/iskip
     write(89,*) nz/iskip
     write(89,*) ny
     write(89,*) xend
     write(89,*) yend
     write(89,*) zend
     write(89,*) c1y
     write(89,*) xloc
     write(89,*) xend
     close(89)

! Write out x & y physical coordinates
!    open (UNIT=93,FILE=cfile,STATUS='UNKNOWN')
!    do j=drange(5),drange(6)
!       do k=drange(3),drange(4),iskip
!          do i=drange(1),drange(2),iskip
!             xx = x(i)
!             zz = z(k)
!             yy = y(j) + phi(i,k)
!             write(93,900) xx, zz, yy
!          end do
!       end do
!    end do
!    close(93)

! Write out coordinates
!    open (UNIT=91,FILE=xcfile,STATUS='UNKNOWN')
!    do i = drange(1), drange(2)
!       write(91,901) i, x(i)
!    end do
!    close(91)

!    open (UNIT=90,FILE=zcfile,STATUS='UNKNOWN')
!    do k = drange(3), drange(4)
!       write(90,901) k, z(k)
!    end do
!    close(90)

!    open (UNIT=92,FILE=ycfile,STATUS='UNKNOWN')
!    do j = drange(5), drange(6)
!       write(92,901) j, y(j), detady(j), deta2dy2(j)
!    end do
!    close(92)

! Write out the surface data...
!    open (UNIT=95,FILE=sfile,STATUS='UNKNOWN')
!    do k = drange(3), drange(4)
!       do i = drange(1), drange(2)
!          write(95,901) x(i), z(k), rb(i,k), phi(i,k),& 
!               dphidx(i,k),dphi2dx2(i,k),psi(i,k,0),&
!               f(i,k,0,neqmax),vvel(i,k,0)
!       end do
!    end do
!    close(95)

!! AI 
!!***************************************************
!!***************************************************
!! GETTING REACTION RATES MANUALLY HERE *************
!!***************************************************
!!***************************************************


     poo = press

     pre(1) = da(1)*poo**np(1)
     pre(2) = da(2)*poo**np(2)
     pre(3) = da(3)*poo**np(3)



     do m=1,3
        do j = drange(5), drange(6)
        do k = drange(3), drange(4)
        do i = drange(1), drange(2)
           rate_out(i,k,j,m) = pre(m)*sign(1.0d0,f(i,k,j,2))*abs(f(i,k,j,2))**mu(m,2)* &
                  sign(1.0d0,f(i,k,j,3))*abs(f(i,k,j,3))**mu(m,3)*sign(1.0d0,f(i,k,j,2))&
                  *abs(f(i,k,j,4))**mu(m,4)*EXP(-thetay(m)/f(i,k,j,1))
        end do
        end do
        end do
     end do


! Write out field variable data...
     open (UNIT=10,FILE=fieldsfile,STATUS='UNKNOWN')
     do j = drange(5), drange(6)
        do k = drange(3), drange(4), iskip
           do i = drange(1), drange(2), iskip
              xx = x(i)
              zz = z(k)
              yy = y(j) + phi(i,k)
              rmax = rate(i,k,j,1)*1.0d-6
              write(10,900) xx, zz, yy, &
                   (f(i,k,j,eqn), eqn=1,neqmax),&
                   rmax, vvel(i,k,j),&
                   phi(i,k), psi(i,k,j), rb(i,k), tout,&
                   rate_out(i,k,j,1),rate_out(i,k,j,2),rate_out(i,k,j,3)
           end do
        end do
     end do
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
  use implicit
  Use Restart_dump

  IMPLICIT NONE

! Local variables
  integer :: inunit,i
  CHARACTER :: comt*100
  REAL*8 ::  ru,lewis_y
!-------------------------------


  nyprocs = 1

  !open(inunit,FILE='input.txt',STATUS='old')
  inunit = 120
  open(inunit,file=inputfile_name)

! Read main inputs
! General Code Parameters
  read(inunit,*) comt
  read(inunit,*) comt
  read(inunit,*) comt
  read(inunit,*) irestart
  read(inunit,*) restartfile
  read(inunit,*) outputfile
  read(inunit,*) npts
  read(inunit,*) timestep
  read(inunit,*) tstop
  read(inunit,*) nyloc,inyloc
  read(inunit,*) writemod
  read(inunit,*) ndump,nsfxrstrt
  read(inunit,*) frame
  read(inunit,*) ncyc_debug,Tncyc,Tjumps,Tfrac
  read(inunit,*) ipack
  read(inunit,*) ncyc_oseen
  read(inunit,*) nperiod
  read(inunit,*) period_begin
  read(inunit,*) ncyc_steady
  read(inunit,*) factorRB
  read(inunit,*) ncells
  read(inunit,*,err=99) sfx,sfxrestart
  goto 101
99 continue
  sfx = 'dat'
  sfxrestart = '.ini'
101 continue
  if(sfxrestart == '.INI') then
     sfxtop = '.TOP'
  else
     sfxtop = '.top'
  endif

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
     STOP 'ERROR OLD INPUT DECK'
  endif

! Equivalence ratio and Homginization parameter
  read(inunit,*) comt
  read(inunit,*) equivalence_ratio
  read(inunit,*) alpha
  read(inunit,*) iwflag



! Solid Phase Physical Parameters
  read(inunit,*) comt
  read(inunit,*) tcold

! Time delays
  read(inunit,*) comt
  read(inunit,*) time_steady_pressure
  read(inunit,*) time_oseen

  !>>> vic: shear, from juPack
  ! Shear 
  READ (inunit, *) comt
  READ (inunit, *) shearFreq
  READ (inunit, *) shearRateAmp
  if (myid == 0) then
     write(6,*)'ShearFreq=',shearFreq
     write(6,*)'ShearRateAmplitude=',shearRateAmp
  endif

!KrierFACTOR
  read(inunit,*)comt
  read(inunit,*) KrierFACTOR
  read(inunit,*) scase
  read(inunit,*) PRESSVARY
  read(inunit,*) HEATVARY
  read(inunit,*) NORADIATION
  read(inunit,*) NOQ4

  read(inunit,*) ntrck
  if(ntrck > 0) then
     allocate(mtrck(ntrck),itrck(ntrck),ktrck(ntrck),jtrck(ntrck))
     do i = 1,ntrck
        read(inunit,*) mtrck(i),itrck(i),ktrck(i),jtrck(i)
     end do
  end if
  nprobes = 0d0


! Close file
  close(inunit)

  call choose_restart_file


  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!............................................................................
! CALCULATE RELEVANT PARAMETERS
!............................................................................
  ru = 1.9859d0
  cp = 0.3d0
  pr = 0.72d0        !Prandtl number

  dim_fact = ru/26.0d0 * 4186.8d0/101325.0d0*1000.0d0 ! kcal/m^3->atm
  dim_fact=1.0d0/dim_fact  !Dimensional factor

  rhoREF = dim_fact * press*1d-3

  pi = ACOS(-1.0d0)
  epi_p=0d0
  omg_p = 0d0
  period_p=1d0
  period_p = max(period_p,1d-12)

! Convert microns to cm
  xend = xend/1.0d4
  yend = yend/1.0d4
  zend = xend

!    Set-up the Lewis numbers
  lewis=1.0d0;

!
!this is the accelration of gravity which should be set to zero if not includede
!
  g_accel = zero


  frameCUT = frame

  open(inunit,FILE='implicit.txt',STATUS='old')
  read(inunit,*) TOLERANCEGAUSS 
  read(inunit,*) iterLINEARsys 
  read(inunit,*) iterativemethod 
  read(inunit,*) GAUSS_GMRES 
  read(inunit,*) implicitBC 
  read(inunit,*) GMRES_restarts 
  close(inunit)
  iterativemethod = 1     ! use gmres, only for radiation
  if (NORADIATION) iterativemethod = 2  ! use newton otherwise

  if(myid == 0) then
     write(*,'(80("&")/,a,/,20("-"))')'output from  INPUTS:'  
     write(*,*) 'KrierFACTOR ',KrierFACTOR
     write(*,*) 'scase ',scase
     write(*,*) 'PRESSVARY ',PRESSVARY
     write(*,*) 'HEATVARY ',HEATVARY
     write(*,*) 'NORADIATION ',NORADIATION
     write(*,*) 'NOQ4 ',NOQ4
     write(*,*) 'Tracking # points',ntrck
     do i = 1,ntrck
        write(*,*) mtrck(i),itrck(i),ktrck(i),jtrck(i)
     end do
     write(*,'(20("-"),/a/,80("&"))')'end INPUTS'
  endif
  

!----------------------------------------------------
  RETURN
END SUBROUTINE INPUTS
!********************************************************



!********************************************************

SUBROUTINE FLUSHH(funit,fname,ll)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

! Local variables
  integer :: funit,ll
  CHARACTER(LEN=ll) fname
!----------------------------------------


  write(*,*) ncyc,'FLUSHING unit',funit,'FILNAM',fname
  close (funit)

  open(UNIT=funit,FILE=trim(fname),POSITION='append')
  write(*,*) ncyc,'FLUSHed'

!----------------------------------------------------
  RETURN
END SUBROUTINE FLUSHH
