Module RESTART_DUMP

  real*8,POINTER :: yst(:),ynow(:),tmpY1(:),tmpY2(:)
  real*8,Allocatable :: hsp(:),alp(:),csp1(:),csp2(:),csp3(:),vsp1(:),vsp2(:),vsp3(:)
Contains  
! ********************************************************************
  SUBROUTINE RESTART(inpname)

    USE GLOBAL_DATA
    USE MYMPI
    use statistics
    IMPLICIT NONE
    include "parallel.h"

!---------------------------------------------------------------------
! Local Variables
    INTEGER ::  P_mpierr, r_unit, old_prm(4),t_unit,nproc_old
    INTEGER ::  i,k,j
    REAL*8 :: dx_old,dz_old,dy_old
    CHARACTER*80 filnam
    CHARACTER(len=*) :: inpname
    LOGICAL log
!---------------------------------------------------------------------

    P_mpierr = 0
    r_unit = 77
    t_unit = 99

    WRITE(filnam,11)TRIM(inpname),myid,SfxRestart
11  FORMAT(a,'_',i3.3,a)


    if(SameTopology) then
       nproc_old = nproc
       dx_old = dx
       dz_old = dz
    else
       open(t_unit,file=TRIM(inpname) // sfxtop,form='formatted',action="read")
       read(t_unit,*)nproc_old
       read(t_unit,*)dx_old,dz_old,dy_old
    endif

    q = zero
    f= zero
    oldsoln = zero
    phi = zero
    pold = zero

    if(abs(dx_old-dx)<1d-8 .and. abs(dz_old - dz)<1d-8 .and. abs(dy_old - dy)<1d-8) then
       if(nproc_old == nproc) then
          call restartSameProcNum
       else
          call restartDiffProcNum
       endif
    else
       if(myid == 0) write(*,*)'Different grids where used, dx/dx_old=',dx/dx_old,&
            &'dx/dx_old=',dz/dz_old,'dy/dy_old=',dy/dy_old
       call restartDiffGrid
    endif
    if(nproc_old /= nproc) close(t_unit)

    do i = drange(1),drange(2)
       do k = drange(3),drange(4)
          if(all(f(i,k,:,:) == 0d0)) then
             write(6,*)'Myid',myid,'Has zeros',i,k,i+ib_mg(1),k+kb_mg(1)
             write(6,*)'FATAL ERROR STOPPING'
             call mpi_finalize(mpi_comm_world)
             stop
          endif
          if(all(q(i,k,:,1:3)== 0d0)) then
             write(6,*)'Myid',myid,'Has zeros Q',i,k,i+ib_mg(1),k+kb_mg(1)
             write(6,*)'WARNING *******************************************'
             write(6,*)'WARNING *******************************************'
             write(6,*)'WARNING *******************************************'
             write(6,*)'WARNING *******************************************'
             write(6,*)'WARNING *******************************************'
          endif
       enddo
    enddo


    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

    p = zero
    pold = zero
    if(myid ==0)write(6,*)'WARNING ZEROING OUT PRESSURE PERTURB in RESTART'

!IF YOU START A NEW SIMULATION USING A BETTER GUESS USE
!    phi_begin = phi(0,0)
!    tinit = 0.0d0
!    phi_per = phi_begin
!    t_per = tinit
!OR TO CONTINUE
    phi_begin = int(phi(0,0)/period + 0.9d0)*period


    if(myid ==0) write(6,*)ncyc, 'RESTARTING FROM TIME ',&
         tinit, t_per,' sec.','PHI = ', phi_begin, phi_per,&
         dint(phi(0,0)/period) - dint(phi_begin/period),&
         dint(phi_begin/period)

    if(ipack == 1) THEN
       CALL QUADRANTS
       CALL  UPDATE_PSI(zero,zero)
    ELSE
       psi = -1.0d0     
    ENDIF
    CALL SOLID_PROP
    CALL PYROLYSIS

!    phi_begin = phi(0)
    phi_begin = zero

    CALL UPDATE_PHI(zero,zero,0)    !evauate the metric derivatives
    !!  CALL BCQ                      !update convection velocity at the boundary

    CALL PRN_GLOB_VARS

!canc---canc=----canc
    !!$  oldvals = 0d0
    !!$  iptroldvals = 0
    !!$  ncyc_stat = 0
    !!$  iptroldvals = 0
    !!$  oldvals = 0d0
    !!$  ncyc_history = 0
    !!$  stat_rb = 0
    !!$  stat_Ts = 0
    !!$  stat_reached = .false.


    RETURN 
!---------------------------------------------------------------------
  Contains


    subroutine restartSameProcNum


      if(myid == 0) write(6,*)'Same Number of Processors'
      INQUIRE (file=TRIM(filnam),exist=log)
      IF (log .EQV. .FALSE.) THEN
         WRITE(6,*)'Solution file for proc ',myid,' does not exist - STOP!'
         WRITE(6,*)'file name: ',filnam
         CALL MPI_Abort (MPI_COMM_WORLD, 1, P_mpierr)
         STOP
      END IF

      OPEN (r_unit,file=TRIM(filnam),form='unformatted',action="read")

      READ(r_unit,err=10,end=10) ncyc, tinit,phi_per,timeStepOld,old_prm,nstat
      call stat_allocate
      READ(r_unit,err=10,end=10) ncyc_history,iptroldvals,stat_rb,stat_Ts,istat,ncyc_stat,stat_reached,oldvals

      log = abs(maxval(old_prm-mg_prm)) > 0
      IF (log .EQV. .TRUE.)THEN
         WRITE(6,*)'OLD parameters differ from new'
         WRITE(6,*)'MG_parameters ',old_prm, mg_prm
         WRITE(6,*)'Continuing.........'
      END IF

      READ(r_unit,err=10,end=10) phi,f,oldsoln(:,:,:,1:maxsize),q!,pold
!!>      READ(r_unit) phi,f,oldsoln(:,:,:,1:maxsize),q!,pold

      CLOSE(r_unit)

      return

10    continue
      write(6,*)myid,'Error in reading the input files(SameProcNum): ',filnam,trim(sfxrestart)
      stop

    end subroutine restartSameProcNum
!*******************************************************


!*******************************************************
    subroutine  restartDiffProcNum
      INTEGER,allocatable :: a_11(:),b_11(:)
      REAL*8, ALLOCATABLE :: phist(:,:),fst(:,:,:,:),qst(:,:,:,:),oldsolnst(:,:,:,:),poldst(:,:,:)
      INTEGER :: i,j,k,ist,kst,id,xrst(3),col(7),myidOld
      include "parallel.h"
!---------------------------------------------------

      allocate(a_11(11))
      allocate(b_11(11))
      b_11(1:11) = alldrange11(1:11,myid+1)

      do id = 0,nproc_old-1

         read(t_unit,*)myidOld
         read(t_unit,*)a_11(1:11)
         read(t_unit,*)filnam
         call check_corners(a_11,b_11,one,one,log)

         if(.not. log) cycle

         write(6,*)'id,myid,myidOld,filnamDiffProcNum,log',id,myid,myidOld,trim(filnam),log
         write(6,*)'myid',myid,a_11,b_11
         write(6,*)''

         OPEN (r_unit,file=TRIM(filnam),form='unformatted',action="read")

         READ(r_unit,err=10,end=20) ncyc, tinit,phi_per,timeStepOld,old_prm,nstat
         call stat_allocate
         READ(r_unit,err=15,end=20) ncyc_history,iptroldvals,stat_rb,stat_Ts,istat,ncyc_stat,stat_reached,oldvals

         xrst(1:3) = a_11(9:11)

         ALLOCATE(phist(-2:xrst(1)+2,-2:xrst(2)+2),STAT=P_mpierr)
         if(P_mpierr >0) write(6,*)'HECK unable to allocate RESTART1',xrst
         ALLOCATE(fst(-2:xrst(1)+2,-2:xrst(2)+2,-2:xrst(3)+2,maxsize),STAT=P_mpierr)
         if(P_mpierr >0) write(6,*)'HECK unable to allocate RESTART2',xrst
         ALLOCATE(qst(-2:xrst(1)+2,-2:xrst(2)+2,-2:xrst(3)+2,ndim),STAT=P_mpierr)
         if(P_mpierr >0) write(6,*)'HECK unable to allocate RESTART3',xrst
         ALLOCATE(oldsolnst(-2:xrst(1)+2,-2:xrst(2)+2,-2:xrst(3)+2,maxsize),STAT=P_mpierr)
         if(P_mpierr >0) write(6,*)'HECK unable to allocate RESTART4',xrst
         ALLOCATE(poldst(-2:xrst(1)+2,-2:xrst(2)+2,-2:xrst(3)+2),STAT=P_mpierr)
         if(P_mpierr >0) write(6,*)'HECK unable to allocate RESTART5',xrst

         READ(r_unit,err=11,end=21) phist,fst,oldsolnst,qst,poldst
         do i = drange(1),drange(2)
            ist = i - a_11(7) + b_11(7)
            do k = drange(3),drange(4)
               kst = k - a_11(8) + b_11(8)
               if(ist >= a_11(1) .and. ist <= a_11(2) .and.&
                    & kst >= a_11(3) .and. kst <= a_11(4)) then
                  q(i,k,:,:) = qst(ist,kst,:,:)
                  f(i,k,:,:) = fst(ist,kst,:,:)
                  oldsoln(i,k,:,1:maxsize) = oldsolnst(ist,kst,:,:)
                  pold(i,k,:) = poldst(ist,kst,:)
                  phi(i,k)   = phist(ist,kst)
               endif
            enddo
         enddo
         CLOSE(r_unit)
         deallocate(phist,fst,qst,oldsolnst,poldst)
      enddo

      CALL MPI_BARRIER(comm3d, P_mpierr) 
      col(1:7) = (/1,2,1,ndim,drange(5:6),0/) 
      finest_mesh%w => q
      call parallel_swap4(col,finest_mesh)
      col(1:7) = (/1,2,1,maxsize,drange(5),drange(6)+1,1/) 
      finest_mesh%w => f
      call parallel_swap4(col,finest_mesh)
      finest_mesh%w => oldsoln
      call parallel_swap4(col,finest_mesh)
      col(1:7) = (/drange(1:4),0,ny,1/) 
      call parallel_swap(pold,col(1:6),finest_mesh)
      CALL FILL_PHI_GHOST
      CALL MPI_BARRIER(comm3d, P_mpierr) 

      return

10    continue

      write(6,*)myid,'Error in reading the input files(DiffProcNum) Prelims : ',filnam,trim(sfxrestart)
      stop

15    continue
      
      write(6,*)myid,'Error in reading the input files(DiffProcNum) Prelims#2 : ',filnam,trim(sfxrestart)
      call mpi_finalize(ierr)
      stop

20    continue

      write(6,*)myid,'End in reading the input files(DiffProcNum) Prelims : ',filnam,trim(sfxrestart)
      stop
11    continue
      write(6,*)myid,'Error in reading the input files(DiffProcNum) Fields  : ',filnam,trim(sfxrestart)
      stop
21    continue
      write(6,*)myid,'End in reading the input files(DiffProcNum) Fields  : ',filnam,trim(sfxrestart)
      stop


    end subroutine restartDiffProcNum
!******************************************************* 

!*******************************************************
    subroutine  restartDiffGrid
      INTEGER,allocatable :: a_11(:),b_11(:)
      REAL*8, ALLOCATABLE :: phist(:,:),fst(:,:,:,:),qst(:,:,:,:)
      REAL*8, ALLOCATABLE :: oldsolnst(:,:,:,:),poldst(:,:,:)
      INTEGER :: i,j,k,ist(2),kst(2),id,xrst(3),col(7),myidOld
      INTEGER :: imedOld,kmedOld,eqn
      REAL*8 :: epsloc,aa(4),Rist,Rkst,RistM,RkstM,ysign
      include "parallel.h"
!---------------------------------------------------


      if(myid == nproc-1) then
         write(6,*) 'DIFFERENT DX',dx_old,dx,'Myid is',myid
         write(6,*) '-------------'
      endif

      epsloc = 1d-12

      allocate(a_11(11))
      allocate(b_11(11))
      b_11(1:11) = alldrange11(1:11,myid+1)

      do id = 0,nproc_old-1

         read(t_unit,*)myidOld
         read(t_unit,*)a_11(1:11)
         read(t_unit,*)filnam
         call check_corners(a_11,b_11,dx_old,dx,log)

         if(.not. log) cycle

         imedOld = a_11(1) + anint((dble(a_11(2))-dble(a_11(1)))/two)
         kmedOld = a_11(3) + anint((dble(a_11(4))-dble(a_11(3)))/two)

         if(myid == nproc-1) then
            write(6,*)'id,myid,myidOld,filnamDiffGrid,log',id,myid,myidOld,trim(filnam),log
            write(6,*)'myid',myid,a_11,b_11
            write(6,*)''
         endif

         OPEN (r_unit,file=TRIM(filnam),form='unformatted',action="read")

         xrst(1:3) = a_11(9:11)

         allocate(yst(-1:xrst(3)+1),ynow(-1:xr(3)+1),STAT=P_mpierr)
         allocate(hsp(-1:xrst(3)+1),alp(-1:xrst(3)+1),csp1(-1:xrst(3)+1),csp2(-1:xrst(3)+1))
         allocate(csp3(-1:xrst(3)+1),vsp1(-1:xrst(3)+1),vsp2(-1:xrst(3)+1),vsp3(-1:xrst(3)+1))
         READ(r_unit,err=10,end=10) ncyc, tinit,phi_per,timeStepOld,old_prm,nstat,yst

         call stat_allocate
         READ(r_unit,err=10,end=10) ncyc_history,iptroldvals,stat_rb,stat_Ts,istat,ncyc_stat,stat_reached,oldvals


         ALLOCATE(phist(-2:xrst(1)+2,-2:xrst(2)+2),STAT=P_mpierr)
         if(P_mpierr >0) write(6,*)'HECK unable to allocate RESTART1',xrst(1:3)
         ALLOCATE(fst(-2:xrst(1)+2,-2:xrst(2)+2,-2:xrst(3)+2,maxsize),STAT=P_mpierr)
         if(P_mpierr >0) write(6,*)'HECK unable to allocate RESTART2',xrst(1:3)
         ALLOCATE(qst(-2:xrst(1)+2,-2:xrst(2)+2,-2:xrst(3)+2,ndim),STAT=P_mpierr)
         if(P_mpierr >0) write(6,*)'HECK unable to allocate RESTART3',xrst(1:3)
         ALLOCATE(oldsolnst(-2:xrst(1)+2,-2:xrst(2)+2,-2:xrst(3)+2,maxsize),STAT=P_mpierr)
         if(P_mpierr >0) write(6,*)'HECK unable to allocate RESTART4',xrst(1:3)
         ALLOCATE(poldst(-2:xrst(1)+2,-2:xrst(2)+2,-2:xrst(3)+2),STAT=P_mpierr)
         if(P_mpierr >0) write(6,*)'HECK unable to allocate RESTART5',xrst(1:3)
         allocate(tmpY1(-1:xrst(3)+1),STAT=P_mpierr)
         allocate(tmpY2(-1:xr(3)+1),STAT=P_mpierr)

         READ(r_unit,err=10,end=10) phist,fst,oldsolnst,qst,poldst

         ynow(-1:xr(3)+1) = y(-1:xr(3)+1)
         do i = drange(1),drange(2)
            Rist = dble(i + b_11(7))*dx/dx_old - dble(a_11(7))
            RistM = Rist + sign(epsloc,dble(imedOld) - Rist )
            ist = (/floor(RistM),ceiling(RistM) /)
            do k = drange(3),drange(4)
               Rkst = dble(k + b_11(8))*dz/dz_old - dble(a_11(8))
               RkstM = Rkst + sign(epsloc,dble(kmedOld) - Rkst )
               kst = (/floor(RkstM),ceiling(RkstM) /)
               if(Rist >= a_11(1)-epsloc .and. Rist < a_11(2)+1 .and.&
                    Rkst >= a_11(3)-epsloc .and. Rkst < a_11(4)+1) then
                  aa(1) = (kst(2) - Rkst)*(ist(2) - Rist)
                  aa(2) = (kst(2) - Rkst)*(Rist - ist(1))
                  aa(3) = (Rkst - kst(1))*(ist(2) - Rist)
                  aa(4) = (Rkst - kst(1))*(Rist - ist(1))

                  do eqn = lbound(q,4),ubound(q,4)
                     do j = -1,xrst(3)+1
                        tmpY1(j) = aa(1)*qst(ist(1),kst(1),j,eqn)+&
                             &       aa(2)*qst(ist(2),kst(1),j,eqn)+&
                             &       aa(3)*qst(ist(1),kst(2),j,eqn)+&
                             &       aa(4)*qst(ist(2),kst(2),j,eqn)
                     enddo
                     call spline(tmpY1,tmpY2,yst,ynow)

                     do j = -1,xr(3)+1
                        q(i,k,j,eqn) = tmpY2(j)
                     end do
                  end do
                  do eqn = lbound(f,4),ubound(f,4)
                     ysign = merge(1d0,-1d0,eqn<neqgas)
                     yst  = yst*ysign
                     ynow = ynow*ysign
                     do j = -1,xrst(3)+1
                        tmpY1(j) = aa(1)*fst(ist(1),kst(1),j,eqn)+&
                             &       aa(2)*fst(ist(2),kst(1),j,eqn)+&
                             &       aa(3)*fst(ist(1),kst(2),j,eqn)+&
                             &       aa(4)*fst(ist(2),kst(2),j,eqn)
                     enddo
                     call spline(tmpY1,tmpY2,yst,ynow)
                     !!$                     if(i*k > 300) then
                     !!$                        do j = 0,xrst(3)
                     !!$                           if(myid == 0) write(121-eqn,*)yst(j),tmpY1(j)
                     !!$                           if(myid == 0) write(122+eqn,*)ynow(j),tmpY2(j)
                     !!$                        enddo
                     !!$                        close(121-eqn)
                     !!$                        close(122+eqn)
                     !!$                        call mpi_barrier(MPI_comm_world,P_mpierr)
                     !!$                        if(eqn == ubound(f,4))call mpi_finalize(MPI_comm_world,P_mpierr)
                     !!$                        if(eqn == ubound(f,4))stop
                     !!$                     endif
                     do j = -1,xr(3)+1
                        f(i,k,j,eqn) = tmpY2(j)
                     end do

                     do j = -1,xrst(3)+1
                        tmpY1(j) = aa(1)*oldsolnst(ist(1),kst(1),j,eqn)+&
                             &       aa(2)*oldsolnst(ist(2),kst(1),j,eqn)+&
                             &       aa(3)*oldsolnst(ist(1),kst(2),j,eqn)+&
                             &       aa(4)*oldsolnst(ist(2),kst(2),j,eqn)
                     enddo
                     call spline(tmpY1,tmpY2,yst,ynow)
                     do j = -1,xr(3)+1
                        oldsoln(i,k,j,eqn) = tmpY2(j)
                     end do
                     yst  = yst/ysign
                     ynow = ynow/ysign
                  end do
                  do j = -1,xrst(3)+1
                     tmpY1(j) = aa(1)*poldst(ist(1),kst(1),j)+&
                          &       aa(2)*poldst(ist(2),kst(1),j)+&
                          &       aa(3)*poldst(ist(1),kst(2),j)+&
                          &       aa(4)*poldst(ist(2),kst(2),j)
                  enddo
                  call spline(tmpY1,tmpY2,yst,ynow)
                  do j = -1,xr(3)+1
                     pold(i,k,j) = tmpY2(j)
                  end do

                  phi(i,k) = aa(1)*phist(ist(1),kst(1))+&
                       &     aa(2)*phist(ist(2),kst(1))+&
                       &     aa(3)*phist(ist(1),kst(2))+&
                       &     aa(4)*phist(ist(2),kst(2))
               endif
            enddo
         enddo
         CLOSE(r_unit)
         deallocate(yst,ynow,hsp,alp)
         deallocate(csp1,csp2,csp3,vsp1)
         deallocate(vsp2,vsp3,tmpY1,tmpY2)
         deallocate(phist,fst,qst,oldsolnst,poldst)
      enddo

      CALL MPI_BARRIER(comm3d, P_mpierr) 
      col(1:7) = (/1,2,1,ndim,drange(5:6),0/) 
      finest_mesh%w => q
      call parallel_swap4(col,finest_mesh)
      col(1:7) = (/1,2,1,maxsize,drange(5),drange(6)+1,1/) 
      finest_mesh%w => f
      call parallel_swap4(col,finest_mesh)
      finest_mesh%w => oldsoln
      call parallel_swap4(col,finest_mesh)
      col(1:7) = (/drange(1:4),0,ny,1/) 
      call parallel_swap(pold,col(1:6),finest_mesh)
      CALL FILL_PHI_GHOST
      CALL MPI_BARRIER(comm3d, P_mpierr) 

      return

10    continue

      write(*,*)myid,'Error in reading the input files--Diffgrid ',filnam
      stop


    end subroutine restartDiffGrid

  END SUBROUTINE RESTART
! ********************************************************************

! ********************************************************************
  SUBROUTINE choose_restart_file
    USE GLOBAL_DATA
    USE MYMPI
    use implicit
    IMPLICIT NONE
    logical :: logic1,logic2
    integer :: inunit,i
    real*8 :: tag1,tag2
    CHARACTER(len=60) :: fname1,fname2
!-------------------------------

    SameTopology = .false.
    fname1 = 'rocfire' // sfxtop
    fname2 = 'rocfout' // sfxtop

    inunit = 121

    inquire(file = trim(fname1), exist = logic1)
    inquire(file = trim(fname2), exist = logic2)
    

    if(logic1) then
       open(inunit,file= trim(fname1),form='formatted',action="read")
       read(inunit,*)i,tag1
       close(inunit)
    else
       tag1 = 1d99
    endif
    if(logic2) then
       open(inunit,file= trim(fname2),form='formatted',action="read")
       read(inunit,*)i,tag2
       close(inunit)
    else
       tag2 = 1d99
    endif
    if(myid == 0) write(*,*)'Looking for files::',trim(fname1),trim(fname2)
    if(myid == 0) write(*,*)'Exist flag::',logic1,logic2
    if(myid == 0) write(*,*)'Tags:: ',tag1,tag2
    if(irestart > 0) then
       if(.not. logic1 .and. .not. logic2) then
          write(*,*)'NO initalization file existent, Better have the same topology'
          SameTopology = .TRUE.
          goto 1000
       endif
       if(restartfile(1:1) /= 'r') then
          if(restartfile(1:1) == 'l' .or. restartfile(1:1) == 'L') then
             restartfile = merge('rocfire','rocfout',(tag1>tag2 .or. tag2>1d50) .and. tag1<1d50)
          elseif(restartfile(1:1) == 'e' .or. restartfile(1:1) == 'E') then
             restartfile = merge('rocfire','rocfout',(tag1<tag2 .or. tag2>1d50) .and. tag1<1d50)
          else
             write(*,*)'Urecognized file name',restartfile
             write(*,*)'STOPPING'
             call mpi_finalize(mpi_comm_world)
             stop
          endif
       endif
    endif
1000 continue
    if(outputfile(1:1) /= 'r') then
       if(outputfile(1:1) == 'l' .or. outputfile(1:1) == 'L') then
          outputfile = merge('rocfire','rocfout',(tag1>tag2 .and. tag2 < 1d50) .or. tag1>1d50 )
       elseif(outputfile(1:1) == 'e' .or. outputfile(1:1) == 'E') then
          outputfile = merge('rocfire','rocfout',(tag1<tag2 .and. tag2 < 1d50) .or. tag1>1d50 )
       else
          write(*,*)'Urecognized file name',outputfile
          write(*,*)'STOPPING'
          call mpi_finalize(mpi_comm_world)
          stop
       endif
    endif


    return
  end SUBROUTINE choose_restart_file
!*******************************************************

!********************************************************************
  SUBROUTINE DUMP(t,inpname)

    USE GLOBAL_DATA
    Use statistics
    IMPLICIT NONE

!---------------------------------------------------------------------
! Dummy Variables
    REAL *8 ::  t

! Local Variables
    INTEGER :: P_mpierr, w_unit,t_unit, id,val(8),cycslct
    REAL*8  :: ref(8),tagnum
    CHARACTER(len=*) :: inpname
    CHARACTER*80 filnam,sfxini,sfxtopo
!---------------------------------------------------------------------


    call date_and_time(VALUES = val)
    ref = (/3.6d3*24d0*31d0*12d0,3.6d3*24d0*31d0,3.6d3*24d0,0d0,3.6d3,6d1,1d0,1d-3/)
    tagnum = sum(ref*dble(val))

    if(ipack /= 1) RETURN
    if(myid == 0 ) CLOSE(9)
    if(myid ==0)write(*,*)'WARNING ZEROING OUT PRESSURE PERTURB in DUMP'

    !!$    cycslct = NINT( dble(mod(ncyc_run,nsfxrstrt*ndump))/ dble(ndump))  !disabled becaus eof problems at the end of the run when 
!ncyc_run is not multiple of ndump

    P_mpierr = 0
    w_unit = 88 

    !!$    if(cycslct == 0) then
    sfxini = '.ini'
    sfxtopo = '.top'
    !!$    else
    !!$       sfxini = '.INI'
    !!$       sfxtopo = '.TOP'
    !!$    endif

    WRITE(filnam,11)TRIM(inpname),myid,TRIM(sfxini)
11  FORMAT(a,'_',i3.3,a)   

    write(*,*)myid,'DUMPING OUT SOLUTION at ncyc = ',ncyc,trim(filnam)

    close(w_unit)
    OPEN (w_unit,file=TRIM(filnam),form='unformatted',action="write")
    rewind(w_unit)

    WRITE(w_unit) ncyc, t, phi_per, timestep, mg_prm,nstat,y(-1:xr(3)+1)
    WRITE(w_unit) ncyc_history,iptroldvals,stat_rb,stat_Ts,istat,ncyc_stat,stat_reached,oldvals

    WRITE(w_unit) phi,f,oldsoln(:,:,:,1:maxsize),q,pold

    CLOSE(w_unit)

    call MPI_BARRIER(comm3d,ierr)


    t_unit = 88
    close(t_unit)
    if(myid == 0) then
       write(*,*)myid,'Dumping',TRIM(inpname) // sfxtopo, t_unit
       open(t_unit,file=TRIM(inpname) // sfxtopo,form='formatted',action="write")
       write(t_unit,*)nproc,tagnum,val,'ncyc = ',ncyc
       write(t_unit,*)dx,dz,dy,'increments'
       do id = 0,nproc-1
          WRITE(filnam,11)TRIM(inpname),id,TRIM(sfxini)
          write(t_unit,*)id
          write(t_unit,*)alldrange11(1:11,id+1)
          write(t_unit,*)TRIM(filnam)
       enddo
       close(t_unit)
    endif

    if(myid == 0 )&
         open(UNIT=9,FILE=errfile,POSITION='APPEND')

!---------------------------------------------------------------------
    RETURN 
  END SUBROUTINE DUMP
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

    CALL INITIAL_CONDITIONS_PHI

    if(ipack == 1) THEN
       CALL QUADRANTS
       CALL UPDATE_PSI(zero,zero)
    ELSE
       psi = -1.0d0     
    ENDIF
    CALL SOLID_PROP
    CALL PYROLYSIS   


    P_mpierr = 0
    r_unit = 77

    filnam = 'print_1d.dat'

    OPEN (r_unit,file=TRIM(filnam),form='formatted',action="read")

    INQUIRE (file=TRIM(filnam),exist=log)
    IF (log .EQV. .FALSE.) THEN
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
    CALL INITIAL_CONDITIONS_Q

    CALL PRN_GLOB_VARS

!---------------------------------------------------------------------
    RETURN 
  END SUBROUTINE RESTART_1D
!******************************************************************

  subroutine check_corners(inpOld,inpNew,Dold,Dnew,out)

    integer,INTENT(IN) :: inpNew(11),inpOLD(11)
    logical,INTENT(OUT) :: out
    real*8,INTENT(IN) :: Dnew,Dold
    integer:: L
    real*8 :: limA(2,2),Corn(2,4),a_11(11),b_11(11),epsloc
!--------------------------------------
    epsloc = 1d-10
    out = .false.

    do L = 1,2
       if(L == 1) then
          a_11 = dble(inpOld)*Dold
          b_11 = dble(inpNew)*Dnew
       else
          a_11 = dble(inpNew)*Dnew
          b_11 = dble(inpOld)*Dold
       endif
       limA(1:2,1) = (/a_11(7) + a_11(1), a_11(7)+ a_11(2)+1/)
       limA(1:2,2) = (/a_11(8) + a_11(3), a_11(8)+ a_11(4)+1/)

       Corn(1:2,1)  = (/b_11(7)+b_11(1), b_11(8)+b_11(3)/)
       Corn(1:2,2)  = (/b_11(7)+b_11(2)+1, b_11(8)+b_11(3)/)
       Corn(1:2,3)  = (/b_11(7)+b_11(2)+1, b_11(8)+b_11(4)+1/)
       Corn(1:2,4)  = (/b_11(7)+b_11(1), b_11(8)+b_11(4)+1/)


       do i = 1,4
          out = out .or. (product(Corn(1,i)-limA(1:2,1)) <= epsloc .and.  &
               &product(Corn(2,i)-limA(1:2,2)) <= epsloc)
       end do
    enddo

    if(out) return
!check the edges

    out = out .or. (any(abs(Corn(2,1)-limA(1:2,2)) <= epsloc) .and.  &
         (Corn(1,1) - limA(1,1))*(Corn(1,2) - limA(2,1)) <= epsloc)
    out = out .or. (any(abs(Corn(2,4)-limA(1:2,2)) <= epsloc) .and.  &
         (Corn(1,4) - limA(1,1))*(Corn(1,3) - limA(2,1)) <= epsloc)
    out = out .or. (any(abs(Corn(1,2)-limA(1:2,1)) <= epsloc) .and.  &
         (Corn(2,2) - limA(1,2))*(Corn(2,3) - limA(2,2)) <= epsloc)
    out = out .or. (any(abs(Corn(1,1)-limA(1:2,1)) <= epsloc) .and.  &
         (Corn(2,1) - limA(1,2))*(Corn(2,4) - limA(2,2)) <= epsloc)


    if(out) return
!check that it is all contained

    out = out .or. ((Corn(1,1) - limA(1,1))*(Corn(1,2) - limA(2,1)) <= epsloc&
         &   .and.  (Corn(2,2) - limA(1,2))*(Corn(2,3) - limA(2,2)) <= epsloc )

    return
  end subroutine check_corners

  subroutine spline(pp1,pp2,cc1,cc2)

    Use Global_data
    Use MyMpi
    implicit none

    real*8,POINTER:: pp1(:),pp2(:),cc1(:),cc2(:)
    real*8 :: prod,yo,der,dcc,difp,difm
    integer :: i,ilow,iup,j,ixtrap
!-----

    ilow = lbound(cc1,1)+1
    iup =  ubound(cc1,1)-1



    prod = 1.0d0

    do i =ilow,iup-1
       hsp(i) = cc1(i+1) - cc1(i)
    enddo
    alp(0) = 0d0
    do i =ilow+1,iup-1
       j = i-1
       alp(i) = 3.0/hsp(i)*(pp1(i+1) -pp1(i)) - 3.0/hsp(j)*(pp1(j+1) -pp1(j))
    enddo
    vsp1(ilow) = 1.d0
    vsp2(ilow) = 0.d0
    vsp3(ilow) = 0.d0
    do i = ilow+1,iup-1
       vsp1(i) = 2.0d0*(cc1(i+1) -cc1(i-1)) - hsp(i-1)*vsp2(i-1)
       vsp2(i) = hsp(i)/vsp1(i)
       vsp3(i) = ( alp(i) - hsp(i-1)*vsp3(i-1) ) / vsp1(i)
    enddo

    vsp1(iup) = 1.d0
    vsp3(iup) = 0.d0
    csp2(iup) = 0.d0
    do i = iup-1,ilow,-1
       csp2(i) = vsp3(i) - vsp2(i)* csp2(i+1)
       csp1(i) = (pp1(i+1) -pp1(i))/hsp(i) -  hsp(i) * ( csp2(i+1)+2.0*csp2(i) ) / 3.0d0
       csp3(i) = (csp2(i+1)-csp2(i))/(3.0d0*hsp(i))
    enddo

!
!
! 
    do j = 0,ubound(cc2,1)
       i=ilow-1
       ixtrap = 0
       prod = 1.d0
       do while(prod .gt. 0)
          i=i+1
          difm=cc2(j)-cc1(i)
          difp=cc2(j)-cc1(i+1)
          prod = difm*difp
          if(i .ge. iup  .AND. prod .gt. 0)then
!>            write(*,101)'FAILED SEARCH',cc2(j),(cc1(m),m=1,n)
             ixtrap = 1
             prod = -1.0
          endif
       enddo
       if(ixtrap == 0) then
          pp2(j) =  pp1(i) + Csp1(i)*(cc2(j)-cc1(i)) +  Csp2(i)*(cc2(j)-cc1(i))**2 + Csp3(i)*(cc2(j)-cc1(i))**3
       elseif(ixtrap .eq. 1) then
          if(abs(cc2(j)-cc1(iup)) .lt. abs(cc2(j)-cc1(ilow)) ) then
             yo = pp1(iup)
             der = (pp1(iup)-pp1(iup-1))/hsp(iup-1)
             dcc = cc2(j)-cc1(iup)
          else
             yo = pp1(ilow)
             der = (pp1(ilow+1)-pp1(ilow))/hsp(ilow)
             dcc = cc2(j)-cc1(ilow)    
          endif
          pp2(j) = yo+der*dcc
       endif
    end do

    return
  end subroutine spline

!
!
!
  Subroutine alldone

    use mympi
    USE GLOBAL_DATA
    integer :: nfile,inunit
    character(LEN=90) :: fname
    REAL*8  :: iteration_time,tmpTime
    logical :: lgex
!-----------------------------------
    if(mod(ncyc_run-1,1000) == 0) then
       inquire(file = 'maxtime.txt',EXIST=lgex)
       if(lgex) then
          inunit = 122
          open(inunit,FILE='maxtime.txt',STATUS='old')
          read(inunit,*) maxtime
          maxtime = maxtime*6d1**2  !seconds
          close(inunit)
       else
          maxtime = 1d99
       endif
    end if
    if(maxtime>= 1d98)return
    tmpTime = mpi_Wtime() - done_time
    call mpi_allreduce(tmpTime,iteration_time,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
    tmpTime = mpi_Wtime() - startTime
    call mpi_allreduce(tmpTime,elapsedTime,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
    if(myid == 0 .and. mod(ncyc,writemod) == 0) write(*,*) ncyc,'alldone---Times',elapsedTime,maxtime,-3d0*iteration_time,elapsedTime-maxtime > -5d0*iteration_time
    if(elapsedTime-maxtime > -3d0*iteration_time) then
       CALL DUMP(tcyc(2),outputfile)
       open(111,file = 'LatestRestart.dat',form = 'formatted')
       write(111,*) ncyc,tcyc
       close(111)
       write(*,*) 'ALL DONE: ',elapsedTime,maxtime
       call mpi_finalize(ierr)
       stop 'alldone'
    end if
    done_time = mpi_Wtime()

    return
  end Subroutine alldone

end Module RESTART_DUMP
