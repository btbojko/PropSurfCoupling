!-------------------------------------------------------------------------------

MODULE MYMPI

  SAVE

  INCLUDE "mpif.h"
  
  integer :: Prid

  !-------------------------------------------------------------------------------
END MODULE MYMPI


Module Timing
  SAVE!
  real*8 :: timeS(100),Ttime,timeT(100)
  character*20 :: timeSTR(100)
contains
  subroutine takeTime(N)
    Use MYMPI
    integer,Intent(IN) :: N
    integer:: M,l,MPI_ERR
    M = abs(n)
    l = sign(1,n)

    if(m > ubound(TimeS,1)) then
       timeS = 0d0
       timeT = 0d0

    else

       if(l < 0) then
          call mpi_barrier(MPI_COMM_WORLD,MPI_ERR)
          timeT(m) = MPI_WTIME()
       else
          call mpi_barrier(MPI_COMM_WORLD,MPI_ERR)
          timeS(m) =  timeS(m) + MPI_WTIME()-timeT(m)
       endif

    endif

    return
  end subroutine takeTime

  subroutine TimeReport
    Use MYMPI
    integer :: iall,i

    if(Prid /= 0) return
    iall = ubound(timeS,1)
    do i = 1,iall
       if(abs(timeS(i)) < 1d-12) cycle
       write(*,'(i3,2x,a,f12.6)')i,'TIMES',timeS(i)
    end do

    return
  end subroutine TimeReport

end Module Timing
