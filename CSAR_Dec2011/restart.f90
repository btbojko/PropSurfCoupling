! ********************************************************************
SUBROUTINE RESTART(inpname)

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!---------------------------------------------------------------------
! Local Variables
  INTEGER ::  P_mpierr, r_unit, old_prm(4)
  REAL*8 :: tstep_old
  CHARACTER*80 filnam
  CHARACTER(len=*) :: inpname
  LOGICAL log
!---------------------------------------------------------------------

   P_mpierr = 0
   r_unit = 77

   WRITE(filnam,11)TRIM(inpname),myid
11 FORMAT(a,'_',i3.3,'.ini')
   
   
   INQUIRE (file=TRIM(filnam),exist=log)
   IF (log .EQV. .FALSE.) THEN
         WRITE(*,*)'Solution file for proc ',myid,' does not exist - STOP!'
         WRITE(*,*)'file name: ',filnam
         CALL MPI_Abort (MPI_COMM_WORLD, 1, P_mpierr)
         STOP
    END IF

    OPEN (r_unit,file=TRIM(filnam),form='unformatted',action="read")

    READ(r_unit) ncyc, tinit,phi_per,tstep_old,old_prm
    
    log = abs(maxval(old_prm-mg_prm)) > 0
    IF (log .EQV. .TRUE.)THEN
         WRITE(*,*)'OLD parameters differ from new - STOP!'
         WRITE(*,*)'MG_parameters ',old_prm, mg_prm
         CALL MPI_Abort (MPI_COMM_WORLD, 1, P_mpierr)
         STOP
    END IF

    READ(r_unit) phi,f,oldsoln,q,pold

    CLOSE(r_unit)

    write(*,*)myid, 'RESTARTING FROM TIME ', tinit,&
         ' LAST PERIOD INFO ',ncyc,tstep_old,timestep
    if(myid ==0)write(*,*)'WARNING ZEROING OUT PRESSURE PERTURB in RESTART'
    p = 0


    tcyc = tinit

    phi_begin = phi(0,0)

    CALL UPDATE_PHI(zero,zero)    !these two depnd on one another
    CALL UPDATE_PSI(zero,zero)    !reinitialize psi counter, ncells might differ
    CALL LAMBDA(0)                ! to handle the boundary conditions
    CALL UPDATE_PHI(zero,zero)    !evauate the metric derivatives
    !!  CALL BCQ                      !update convection velocity at the boundary

!---------------------------------------------------------------------
  RETURN 
END SUBROUTINE RESTART
!*************************************************************

!********************************************************************
SUBROUTINE DUMP(t,inpname)

  USE GLOBAL_DATA
  IMPLICIT NONE

!---------------------------------------------------------------------
! Dummy Variables
  REAL *8 ::  t

! Local Variables
  INTEGER ::  P_mpierr, w_unit
  CHARACTER(len=*) :: inpname
  CHARACTER*80 filnam
!---------------------------------------------------------------------

  if(myid ==0)write(*,*)'WARNING ZEROING OUT PRESSURE PERTURB in DUMP'
  p = 0

  if(myid == 0 ) CLOSE(9)

  P_mpierr = 0
  w_unit = 88

  WRITE(filnam,11)TRIM(inpname),myid
11 FORMAT(a,'_',i3.3,'.ini')   

   write(*,*)myid,'DUMPING OUT SOLUTION at ncyc = ',ncyc,TRIM(filnam)
   
   OPEN (w_unit,file=TRIM(filnam),form='unformatted',action="write")

   WRITE(w_unit) ncyc, t, phi_per, timestep, mg_prm

   WRITE(w_unit) phi,f,oldsoln,q,pold

   CLOSE(w_unit)

   if(myid == 0 )&
   open(UNIT=9,FILE=errfile,POSITION='APPEND')

!---------------------------------------------------------------------
  RETURN 
END SUBROUTINE DUMP
!*************************************************************


!********************************************************************
SUBROUTINE INIT_FROM_SCRATCH

USE GLOBAL_DATA
IMPLICIT NONE

!---------------------------------------------------------------------
INTEGER :: i,k,j
REAL *8 :: burn_rate, coe,poo,qsur,pressure,rhogas
!---------------------------------------------------------------------

ncyc = 0
tinit = 0.0d0
phi_begin = 0.0d0
t_per = tinit
phi_per = phi_begin

CALL INITIAL_CONDITIONS1

CALL UPDATE_PSI(zero,zero)

CALL INITIAL_CONDITIONS2

if(myid == 0) write(*,*)'STARTING FROM SCRATCH'

!.....................................................................................
!  Compute initial burn rates AND update velocity Boundary Conditions
!.....................................................................................
   CALL PYROLYSIS

   CALL BCQ

!---------------------------------------------------------------------
  RETURN 
END SUBROUTINE INIT_FROM_SCRATCH
!*************************************************************

! ********************************************************************
SUBROUTINE RESTART_1D

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!---------------------------------------------------------------------
! Local Variables
  INTEGER ::  r_unit,P_mpierr, i,k,j
  CHARACTER*80 filnam
  REAL*8  :: dummy
  REAL*8,ALLOCATABLE :: sol1d(:,:)
  LOGICAL log
!---------------------------------------------------------------------

   ALLOCATE(sol1d(0:ny,5))

   tinit = 0.0d0
   phi_begin = 0.0d0
   t_per = tinit
   phi_per = phi_begin
   ncyc = 0

   CALL INITIAL_CONDITIONS1
   CALL UPDATE_PSI(zero,zero)

 
   P_mpierr = 0
   r_unit = 77

   filnam = 'print_1d.dat'

   OPEN (r_unit,file=TRIM(filnam),form='formatted',action="read")
   
   INQUIRE (file=TRIM(filnam),exist=log)
   IF (log .EQV. .FALSE.) THEN
         WRITE(*,*)'Solution file for proc ',myid,' does not exist - STOP!'
         WRITE(*,*)'file name: ',filnam
         CALL MPI_Abort (MPI_COMM_WORLD, 1, P_mpierr)
         STOP
   END IF
    
   DO j = 0,ny
      READ(r_unit,*) dummy,sol1d(j,1:5)
      if(myid == 0)&
           write(*,*)j,'INITIAL SOLUTION',dummy,sol1d(j,1:5)
   ENDDO

   CLOSE(r_unit)

   do k = drange(3)-1, drange(4)+1
      do i = drange(1)-1, drange(2)+1
         do j = 0,ny
            f(i,k,j,1:5) = sol1d(j,1:5) 
         enddo
      end do
   end do
  
   DEALLOCATE(sol1d)

   CALL PYROLYSIS
   if(myid == 0) write(*,*)'STARTING FROM 1 D SOLUTION'
!---------------------------------------------------------------------
  RETURN 
END SUBROUTINE RESTART_1D
!******************************************************************
