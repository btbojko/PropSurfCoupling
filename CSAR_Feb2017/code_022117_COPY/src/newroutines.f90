
!**********************************************************
SUBROUTINE NAN_VAL(ff,outcome,neq)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!-----------------------------------------------------------
!!!   Local Variables
  integer i, k, j, n,neq
  real *8,INTENT(IN) :: ff(neq)
  REAL*8,allocatable ::  mx(:),mn(:)
  LOGICAL :: outcome
!-----------------------------------------------------------

!
  allocate(mx(neq),mn(neq))
  
  mx = 1d5
  mn = -1d5

  outcome = .FALSE.
  do n=1,neq
     outcome = outcome .OR. (.NOT. ff(n) > mn(n)) .OR. (.NOT. ff(n) < mx(n))
  enddo

  deallocate(mx,mn)

!--------------------------------------------------------
  RETURN
END SUBROUTINE NAN_VAL
!**********************************************************************
!

!**********************************************************
SUBROUTINE NAN_FIELD(outcome,maxprint)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!-----------------------------------------------------------
!!!   Local Variables
  INTEGER :: i,k,j,nprinted,maxprint
  LOGICAL :: outcome,repeat,doprint
!-----------------------------------------------------------
!
!!check for NANs
  repeat = .FALSE.
  nprinted = 0;
  do i = drange(1),drange(2)
     do k = drange(3),drange(4)
        do j = drange(5)+1, drange(6)-1

           call NAN_VAL(f(i,k,j,1:neqgas),outcome,neqgas)
           repeat = repeat .or. outcome
          
           if(outcome .and. nprinted <= maxprint)then
              nprinted = nprinted+1
              write(*,'(i8,4i5,a,1p5e12.4,a,5e12.4,a,10e12.4,a,4e12.4)')&
                    ncyc,myid,i,k,j,'OUT OF RNGE F',&
                   f(i,k,j,1:neqgas+1),'OLD',&
                   oldsoln(i,k,j,1:neqgas+1)
           endif
 
           call NAN_VAL(q(i,k,j,1:ndim),outcome,ndim)
           repeat = repeat .or. outcome

           if(outcome .and. nprinted <= maxprint)then
              nprinted = nprinted+1
              write(*,'(i8,4i5,a,1p3e12.4,a,3e12.4,a,10e12.4,a,3e12.4)')&
                    ncyc,myid,i,k,j,'OUT OF RNGE Q',&
                   q(i,k,j,1:ndim),'OLD',&
                   oldsoln(i,k,j,maxsize+1:maxsize+ndim)
           endif

        enddo
     enddo
  enddo

  call MPI_ALLREDUCE(repeat,outcome,1,MPI_LOGICAL,MPI_LOR,comm3d,ierr)

!--------------------------------------------------------
  RETURN
END SUBROUTINE NAN_FIELD
!**********************************************************************


!**********************************************************
SUBROUTINE NAN_DFFIELD(outcome,maxprint,msg)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  !-----------------------------------------------------------
!!!   Local Variables
  INTEGER :: i,k,j,nprinted,maxprint
  LOGICAL :: outcome,repeat,doprint
  character*(*) :: msg
  !-----------------------------------------------------------
  !
  !!check for NANs
  repeat = .FALSE.
  nprinted = 0;
  do i = drange(1)-1,drange(2)+1
     do k = drange(3)-1,drange(4)+1
        do j = drange(5)-1, drange(6)+1

           call NAN_VAL(dff(i,k,j,1:neqgas),outcome,neqgas)
           repeat = repeat .or. outcome

           if(outcome .and. nprinted <= maxprint)then
              nprinted = nprinted+1
              write(*,'(i8,4i5,a,a,1p5e12.4,a,5e12.4,a,10e12.4,a,4e12.4)')&
                   ncyc,myid,i,k,j,'OUT OF RNGE DFF',trim(msg),dff(i,k,j,1:neqgas)
           endif
        enddo
     enddo
  enddo

  call MPI_ALLREDUCE(repeat,outcome,1,MPI_LOGICAL,MPI_LOR,comm3d,ierr)

  !--------------------------------------------------------
  RETURN
END SUBROUTINE NAN_DFFIELD
!**********************************************************************

!**********************************************************
SUBROUTINE NAN_STOP(msg,id)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!-----------------------------------------------------------
!!!   Local Variables
  INTEGER :: i,k,j,id
  LOGICAL :: outcome,repeat,doprint
  character*(*)::msg
!-----------------------------------------------------------

  call NAN_FIELD(outcome,10)
  if(outcome) then
    write(*,*)msg,id
    write(*,*)'SHUTTING OFF BECAUSE OF NON_convergence,myid :: ',myid
    call MPI_ABORT(MPI_COMM_WORLD, i,j)
    call MPI_FINALIZE(i)
    stop
  endif

!--------------------------------------------------------
  RETURN
END SUBROUTINE NAN_STOP
!**********************************************************************



!**********************************************************
SUBROUTINE checkIndex(ii,jj,kk,nn)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!-----------------------------------------------------------
!!!   Local Variables
  INTEGER :: i,k,j,ii(1:nn),jj(1:nn),kk(1:nn),n,nn
  LOGICAL :: outcome,repeat,doprint
!-----------------------------------------------------------


  if(any(ii < jj) .or. any(ii > kk)) then
     write(*,*)'OUT of Bounds array','II',ii,'JJ',jj,'KK',kk
     write(*,*)'##################'
  endif
!--------------------------------------------------------
  RETURN
END SUBROUTINE CHECKINDEX
!**********************************************************************
!
! ********************************************************************                                   
SUBROUTINE PROBESETUP

  USE data_types
  USE GLOBAL_DATA

  IMPLICIT NONE

!-------------------------------------------------------------------------                               
!     Local Variables
!------------------------------------------------------------------

  solidProbes(1)%i = floor(dble(drange(2))*half)
  solidProbes(1)%k = floor(dble(drange(4))*half)
  solidProbes(1)%X = x(solidProbes(1)%i)
  solidProbes(1)%Z = z(solidProbes(1)%k)
  solidProbes(1)%PHI = phi(solidProbes(1)%i,solidProbes(1)%k)
  solidProbes(1)%Y = -300d-4 + solidProbes(1)%PHI
  solidProbes(1)%J = floor(.1159757998d0/dy)
  close(101+myid)
  
  return
end SUBROUTINE PROBESETUP



SUBROUTINE PROBEUPDATE

  USE data_types
  USE GLOBAL_DATA

  IMPLICIT NONE

!-------------------------------------------------------------------------                               
!     Local Variables
  integer :: findj,i,k,jloc_P,eqn
  real*8:: yloc_P,signloc
!------------------------------------------------------------------
  
  i = solidProbes(1)%i;k=solidProbes(1)%k
  yloc_P = abs(solidProbes(1)%Y - phi(solidProbes(1)%i,solidProbes(1)%k))
  signloc = sign(1d0,solidProbes(1)%Y - phi(solidProbes(1)%i,solidProbes(1)%k))
  eqn = merge(1,neqmax,signloc>0d0)

  jloc_P = findj(yloc_P,solidProbes(1)%J)
  
  solidProbes(1)%T = f(i,k,jloc_P,eqn)*(y(jloc_P+1)-yloc_P)/(y(jloc_P+1)-y(jloc_P)) - &
       &f(i,k,jloc_P+1,eqn)*(y(jloc_P)-yloc_P)/(y(jloc_P+1)-y(jloc_P))
  write(101+myid,*) ncyc,int(signloc),tcyc(2),solidProbes(1)%X,solidProbes(1)%Z,&
       &solidProbes(1)%Y,solidProbes(1)%T,psi(i,k,jloc_P)

  solidProbes(1)%J = jloc_P
  if(mod(ncyc,1) == 0) then
     close (101+myid)
     open(UNIT=101+myid,POSITION='append')
  endif
  return
end SUBROUTINE PROBEUPDATE



SUBROUTINE PROBERESTART

  USE data_types
  USE GLOBAL_DATA

  IMPLICIT NONE

!------------------------------------------------------------------------- 
!     Local Variables
  integer :: ncyc_R,i_R,eleread
  real*8 :: tm_r,X_R,Z_R,Y_R,T_R,psi_R
!------------------------------------------------------------------                                      

  eleread = 0
  do while (.true.)
     read(101+myid,*,ERR = 101,END=101) ncyc_R,i_R,tm_r,X_R,Z_R,Y_R,T_R,psi_R
     eleread = eleread+1
  enddo

  101 continue

  if (eleread <=1) then
     call PROBESETUP
     return
  endif

  solidProbes(1)%i = floor(dble(drange(2))*half)
  solidProbes(1)%k = floor(dble(drange(4))*half)
  solidProbes(1)%X = x(solidProbes(1)%i)
  solidProbes(1)%Z = z(solidProbes(1)%k)
  solidProbes(1)%PHI = phi(solidProbes(1)%i,solidProbes(1)%k)
  solidProbes(1)%Y = Y_R
  solidProbes(1)%J = floor(.1159757998d0/dy)
  close (101+myid)
  open(UNIT=101+myid,POSITION='append')
  write(101+myid,*) '%%% NEW RESTART DUMP'
  write(*,*) 'DONE RESTART DUMP'
  close (101+myid)
  open(UNIT=101+myid,POSITION='append')
  
  return
end SUBROUTINE PROBERESTART




integer function findJ(yy,jj)

  
  USE data_types
  USE GLOBAL_DATA

  IMPLICIT NONE
  
  integer :: jj,j
  real*8 :: yy,jR,eps_loc


  eps_loc = 1d-6
  jR = dble(jj)

  do while (abs(ggrid(JR)) > 1d-12)
     Jr = Jr - ggrid(JR)*two*eps_loc/(ggrid(JR+eps_loc)-ggrid(JR-eps_loc))
  enddo

  
  findj = floor(Jr)
 
  
  
  return

  contains 
    real*8 function ggrid(JJ)
      real*8 JJ
      
      ggrid = yy-JJ*dy/((2-JJ**2*dy**2/yend**2)**c1y)

      return
    end function ggrid


  end function findJ
