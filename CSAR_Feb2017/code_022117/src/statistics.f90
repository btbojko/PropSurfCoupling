module statistics
  use global_data
  USE MYMPI

  SAVE

  real*8,allocatable :: stat_rb(:),stat_Ts(:)
  integer,allocatable :: istat(:),ncyc_stat(:)
  logical,allocatable :: stat_reached(:)
  real*4,allocatable :: oldvals(:)
  integer :: iadd = 0,nstat = 7,ncyc_history = 0
  integer :: stat_allocated = 0,iptroldvals
  real*8 :: phiAVG

contains

  subroutine stat_allocate

    if(stat_allocated > 0) return
    
    stat_allocated = 1
    
    nstat = 7
    allocate(stat_rb(nstat),stat_Ts(nstat),istat(nstat),ncyc_stat(nstat),stat_reached(nstat))
    stat_rb = 0
    stat_Ts = 0
    istat = 0
    if(irestart == 0)  ncyc_history = ncyc-1
    stat_reached = .false.
    

    istat = (/ 1, 200, 1000, 5000, 25000, 100000, 1000000/)

    allocate(oldvals(min(maxval(istat),100000)))
    oldvals = 0d0
    iptroldvals = 0
    ncyc_stat = 0
    
  end subroutine stat_allocate

  
  subroutine stat_average

    implicit none

    integer :: i,k,j,m,kk
    real*8 :: tmpval(4)
    real*8 :: allval(4)
    real*8 :: const
    !------------------------------------------------

    !average on the surface

    call stat_allocate
    ncyc_history = ncyc_history + 1

    tmpval = 0d0
    do k=drange(3),drange(4)
       do i=drange(1),drange(2)
          const = sqrt(1.0+dphidx(i,k)**2+dphidz(i,k)**2)
          tmpval = tmpval + (/ phit(i,k), f(i,k,0,1),phi(i,k),1d0 /)
       enddo
    enddo

    call mpi_allreduce(tmpval,allval,4, MPI_Double_Precision, MPI_SUM, MPI_COMM_WORLD,ierr) 

    allval(1:3) = allval(1:3)/allval(4)
    m = 1
    stat_rb(m) = allval(1)
    stat_Ts(m) = allval(2)
    phiAVG = allval(3)
    ncyc_stat(m) = 1
    stat_reached(m) = .true.
    do m = 2,ubound(istat,1)
       stat_rb(m) = stat_rb(m) + allval(1)
       stat_Ts(m) = stat_Ts(m) + allval(2)
       ncyc_stat(m) = ncyc_stat(m) + 1
       if(ncyc_stat(m) > istat(m)) then
          if( abs(oldvals(ncyc_history - istat(m) - iptroldvals)) > 1d-14) then
             stat_rb(m) = stat_rb(m) - oldvals(ncyc_history - istat(m) - iptroldvals)
          else
             stat_rb(m) = stat_rb(m) - oldvals(count(abs(oldvals(1:ncyc_history - iptroldvals)) < 1d-14) +1)
          endif
          ncyc_stat(m) = min(ncyc_stat(m),istat(m))
          stat_reached(m) = .true.
       endif
    enddo
    if(ncyc_history - iptroldvals > ubound(oldvals,1)) then
       iptroldvals = iptroldvals + 1
       oldvals (1) = allval(1)
       oldvals = cshift(oldvals,1)
    else
       oldvals(ncyc_history - iptroldvals) = allval(1)
    endif

  end subroutine stat_average

end module statistics
